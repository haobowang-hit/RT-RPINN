"""
3-Head LSTM PINN for Shape Memory Materials
Example 1: 3D Solid Block - Uniaxial Tension and Shape-Memory Recovery

Modification from lstm_pinn_ex1.py:
- LSTM_PINN_3Head replaces LSTM_PINN.
  Three independent LSTM branches, one per displacement component:
    Head U1: spatial encoder -> LSTM -> decoder -> u1  (main loading direction)
    Head U2: spatial encoder -> LSTM -> decoder -> u2  (coupling direction)
    Head U3: spatial encoder -> LSTM -> decoder -> u3  (thickness / Poisson)
  Internal variables q predicted from the U1 LSTM features (-> 6 x N_prony).

  Motivation: in the shared-LSTM baseline the dominant U1 gradient (~28x larger
  than U3) suppresses the U3 signal during backpropagation.  Separate heads
  remove that cross-component interference.

- FEDataLoader, MaterialParameters, LSTM_PINNSolver (all loss functions and the
  training loop) are imported unchanged from lstm_pinn_ex1.py.
- Console output is mirrored to lstm_pinn_3head_EX1.txt.
"""

import sys
import numpy as np
import torch
import torch.nn as nn
from pathlib import Path

torch.manual_seed(42)
np.random.seed(42)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# ── Shared components (no changes) ───────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))

from spatiotemporal_pinn_ex1 import (
    FEDataLoader,
    MaterialParameters,
    plot_training_history,
    plot_displacement_field,
    plot_shape_memory_cycle,
)
from lstm_pinn_ex1 import LSTM_PINNSolver


# ── Console / file tee ────────────────────────────────────────────────────────
class Tee:
    """Write output to multiple streams simultaneously."""

    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()


# ── 3-Head LSTM Architecture ──────────────────────────────────────────────────
class LSTM_PINN_3Head(nn.Module):
    """
    3-Head LSTM PINN: separate LSTM branch for each displacement component.

    Architecture per head (U1 / U2 / U3) - identical hyper-parameters:
      Input (x, y, z, t, T)
        -> Spatial encoder  (MLP with Tanh, independent weights per head)
        -> LSTM             (temporal dependency, independent per head)
        -> Decoder          (MLP with Tanh, outputs 1 scalar displacement)

    Internal variables q:
      Produced by a linear layer on the U1 LSTM hidden state.
      q shape: (batch, 6, N_prony).  Same interface as the baseline model.

    Parameters
    ----------
    spatial_layers : list
        Layer widths for the spatial encoder, e.g. [5, 256, 256, 256].
        First element must be 5 (x, y, z, t, T).
    lstm_hidden : int
        Number of LSTM hidden units per head.
    lstm_layers : int
        Number of stacked LSTM layers per head.
    output_layers : list
        Decoder hidden widths + output dim.  Last element must be 1, e.g. [256, 128, 1].
    n_prony : int
        Number of Prony series branches.
    output_internal_vars : bool
        If True, the forward pass returns q as a 4th output.
    dropout : float
        Dropout rate applied after each encoder layer and between LSTM layers.
    """

    def __init__(self, spatial_layers=[5, 256, 256, 256], lstm_hidden=512, lstm_layers=3,
                 output_layers=[256, 128, 1], n_prony=6, output_internal_vars=True, dropout=0.02):
        super(LSTM_PINN_3Head, self).__init__()

        self.n_prony              = n_prony
        self.output_internal_vars = output_internal_vars
        self.lstm_hidden          = lstm_hidden
        self.lstm_layers          = lstm_layers

        # ------------------------------------------------------------------
        # Builder helpers
        # ------------------------------------------------------------------
        def _make_encoder(layers, drop):
            modules = []
            for i in range(len(layers) - 1):
                modules.append(nn.Linear(layers[i], layers[i + 1]))
                modules.append(nn.Tanh())
                if drop > 0:
                    modules.append(nn.Dropout(drop))
            return nn.Sequential(*modules)

        def _make_decoder(hidden_size, out_layers):
            dims    = [hidden_size] + out_layers
            modules = []
            for i in range(len(dims) - 1):
                modules.append(nn.Linear(dims[i], dims[i + 1]))
                if i < len(dims) - 2:          # no activation after the final layer
                    modules.append(nn.Tanh())
            return nn.Sequential(*modules)

        lstm_kwargs = dict(
            input_size  = spatial_layers[-1],
            hidden_size = lstm_hidden,
            num_layers  = lstm_layers,
            batch_first = True,
            dropout     = dropout if lstm_layers > 1 else 0,
        )

        # ------------------------------------------------------------------
        # Head U1  (main loading direction)
        # ------------------------------------------------------------------
        self.encoder_u1 = _make_encoder(spatial_layers, dropout)
        self.lstm_u1    = nn.LSTM(**lstm_kwargs)
        self.decoder_u1 = _make_decoder(lstm_hidden, output_layers)

        # ------------------------------------------------------------------
        # Head U2  (coupling direction)
        # ------------------------------------------------------------------
        self.encoder_u2 = _make_encoder(spatial_layers, dropout)
        self.lstm_u2    = nn.LSTM(**lstm_kwargs)
        self.decoder_u2 = _make_decoder(lstm_hidden, output_layers)

        # ------------------------------------------------------------------
        # Head U3  (thickness / Poisson direction)
        # ------------------------------------------------------------------
        self.encoder_u3 = _make_encoder(spatial_layers, dropout)
        self.lstm_u3    = nn.LSTM(**lstm_kwargs)
        self.decoder_u3 = _make_decoder(lstm_hidden, output_layers)

        # ------------------------------------------------------------------
        # Internal variables q  (linear layer on U1 LSTM features)
        # ------------------------------------------------------------------
        if output_internal_vars:
            self.internal_var_decoder = nn.Linear(lstm_hidden, 6 * n_prony)

        self._initialize_weights()

    # ------------------------------------------------------------------
    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)

    # ------------------------------------------------------------------
    def forward(self, x, y, z, t, T, return_hidden=False):
        """
        Parameters
        ----------
        x, y, z, t, T : Tensor  shape (B, 1)
        return_hidden  : bool  – also return U1 LSTM feature vector

        Returns
        -------
        u1, u2, u3 : Tensor  shape (B, 1)
        q          : Tensor  shape (B, 6, N_prony)  [only if output_internal_vars]
        feat1      : Tensor  shape (B, lstm_hidden)  [only if return_hidden]
        """
        inp = torch.cat([x, y, z, t, T], dim=1)   # (B, 5)

        # Head U1
        sf1          = self.encoder_u1(inp)                  # (B, enc_dim)
        out1, _      = self.lstm_u1(sf1.unsqueeze(1))        # (B, 1, lstm_hidden)
        feat1        = out1[:, -1, :]                         # (B, lstm_hidden)
        u1           = self.decoder_u1(feat1)                 # (B, 1)

        # Head U2
        sf2          = self.encoder_u2(inp)
        out2, _      = self.lstm_u2(sf2.unsqueeze(1))
        feat2        = out2[:, -1, :]
        u2           = self.decoder_u2(feat2)

        # Head U3
        sf3          = self.encoder_u3(inp)
        out3, _      = self.lstm_u3(sf3.unsqueeze(1))
        feat3        = out3[:, -1, :]
        u3           = self.decoder_u3(feat3)

        if self.output_internal_vars:
            q = self.internal_var_decoder(feat1).reshape(-1, 6, self.n_prony)
            if return_hidden:
                return u1, u2, u3, q, feat1
            return u1, u2, u3, q
        else:
            if return_hidden:
                return u1, u2, u3, feat1
            return u1, u2, u3

    # ------------------------------------------------------------------
    def count_parameters(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    script_dir = Path(__file__).parent
    log_path   = script_dir / 'lstm_pinn_3head_EX1.txt'
    log_fh     = open(log_path, 'w', encoding='utf-8')
    sys.stdout = Tee(sys.__stdout__, log_fh)

    try:
        print("=" * 70)
        print("3-Head LSTM-PINN  |  Example 1  |  Shape Memory Recovery")
        print("Architecture: 3 independent LSTM branches (U1 / U2 / U3)")
        print("=" * 70)
        print()

        # ── FE data ───────────────────────────────────────────────────────
        print("Loading FE data...")
        print("-" * 70)
        fe_loader = FEDataLoader(
            script_dir / 'EX-1-RESULTS',
            script_dir / 'step-frame-time.csv',
        )
        fe_loader.load_all_data()
        bounds = fe_loader.get_domain_bounds()
        print("\nDomain bounds:")
        for key, val in bounds.items():
            print(f"  {key}: {val:.2f}")
        print()

        # ── Material ──────────────────────────────────────────────────────
        print("Initializing material parameters...")
        print("-" * 70)
        mat = MaterialParameters()
        print(f"  N_prony  : {mat.N_prony}")
        print(f"  C11_inf  : {mat.C11_inf / 1e9:.2f} GPa")
        print(f"  T_ref    : {mat.T_ref} K  ({mat.T_ref - 273.15:.1f} °C)")
        print(f"  T_switch : {mat.T_switch} K  ({mat.T_switch - 273.15:.1f} °C)")
        print(f"  WLF      : C1={mat.C1}, C2={mat.C2}")
        print()

        # ── Model ─────────────────────────────────────────────────────────
        print("Initializing 3-Head LSTM-PINN model...")
        print("-" * 70)
        model = LSTM_PINN_3Head(
            spatial_layers    = [5, 256, 256, 256],
            lstm_hidden       = 512,
            lstm_layers       = 3,
            output_layers     = [256, 128, 1],
            n_prony           = mat.N_prony,
            output_internal_vars = True,
            dropout           = 0.02,
        )
        print("Model architecture (3-Head LSTM-PINN):")
        print("  [Head U1]  Encoder [5,256,256,256]  LSTM 512x3  Decoder [256,128,1]")
        print("  [Head U2]  Encoder [5,256,256,256]  LSTM 512x3  Decoder [256,128,1]")
        print("  [Head U3]  Encoder [5,256,256,256]  LSTM 512x3  Decoder [256,128,1]")
        print(f"  [q-head ]  Linear(512 -> 6x{mat.N_prony}) on U1 LSTM features")
        print(f"  Total parameters: {model.count_parameters():,}")
        print()

        # ── Solver (unchanged) ────────────────────────────────────────────
        solver = LSTM_PINNSolver(model, mat, fe_loader, bounds, use_physics_loss=True)
        print("Solver initialized (LSTM_PINNSolver, unchanged from lstm_pinn_ex1).")
        print(f"  L_ref  = {solver.L_ref} mm   |   E_ref = {solver.E_ref / 1e9:.1f} GPa")
        print(f"  Weights: data={solver.lambda_data}, pde={solver.lambda_pde}, "
              f"bc={solver.lambda_bc}, evol={solver.lambda_evol}")
        print(f"  U-weights: U1={solver.lambda_u1}, U2={solver.lambda_u2}, "
              f"U3={solver.lambda_u3}")
        print()

        # ── Training ──────────────────────────────────────────────────────
        print("=" * 70)
        print("Training")
        print("=" * 70)
        solver.train(
            epochs       = 10000,
            batch_size   = 1024,
            learning_rate= 1e-4,
            n_pde_points = 1000,
            n_bc_points  = 200,
        )
        solver.save_loss_history(script_dir / 'lstm_3head_EX1_loss_history.csv')
        print()

        # ── Visualisations ────────────────────────────────────────────────
        print("Generating visualizations...")
        print("-" * 70)
        plot_training_history(solver, save_dir=script_dir)
        plot_displacement_field(solver, fe_loader, save_dir=script_dir)
        plot_shape_memory_cycle(solver, fe_loader, save_dir=script_dir)
        print()

        # ── Save model ────────────────────────────────────────────────────
        model_path = script_dir / 'lstm_pinn_3head_EX1_model.pth'
        torch.save(model.state_dict(), model_path)
        print(f"Model saved : {model_path}")
        print(f"Log saved   : {log_path}")
        print()

    finally:
        log_fh.close()
        sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()
