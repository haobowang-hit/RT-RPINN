# SMPCPINN: SMPC-Informed Spatiotemporal PINN for Paired Forward–Inverse Tasks

This repository provides a physics-informed learning workflow for **history-dependent thermo-viscoelastic shape memory polymer composites (SMPCs)**. The distinguishing feature of SMPCPINN is that it is designed around a **paired forward–inverse task setting**:

* **Forward task:** predict/reconstruct the spatiotemporal response fields (displacement, and derived strain/stress) under prescribed loading and temperature schedules.
* **Inverse task:** identify constitutive parameters from sparse observations (field samples and/or boundary observables such as reaction force), while enforcing the same governing constraints used in the forward task.

The benchmark suite contains **four representative examples ** that form a progressive validation ladder from homogeneous to heterogeneous fields, and from isothermal to non-isothermal identification.

---

## What SMPCPINN Solves

SMPC behaviors are strongly **history-dependent** and often governed by **time–temperature shifting**. High-fidelity FE simulations can resolve these effects but become expensive when the time horizon is long, the temperature schedule is multi-step, or repeated evaluations are needed for calibration and design. SMPCPINN targets a practical balance:

* Data efficiency: only sparse reference supervision is required.
* Physics consistency: governing constraints and boundary conditions are enforced during training.
* Dual capability: the same framework supports both forward prediction and inverse identification.

---

## Two Base Models: PINNBase vs LSTMBase

SMPCPINN provides two interchangeable temporal backbones. They share the same training interface and loss structure; only the temporal representation differs.

### PINNBase (MLP-based spatiotemporal PINN)

* Uses a conventional spatiotemporal neural field.
* Treats time (and temperature) as regular inputs.
* Strong baseline for many forward problems and shorter/less stiff histories.

### LSTMBase (sequence-aware PINN)

* Replaces the temporal representation with an LSTM-based sequence encoder.
* Explicitly propagates temporal state across frames, improving stability for long-range memory effects.
* Particularly beneficial when:

  * the response contains long relaxation tails,
  * temperature schedules induce stiff time–temperature mapping,
  * inverse identification relies on stable strain increments and consistent history propagation.

In the benchmark suite, Ex1–Ex4 can be run with either PINNBase or LSTMBase to quantify the impact of sequence-aware encoding.

---

## Benchmark Suite : Tasks and Design Logic

The four examples are intentionally designed to cover:

1. forward evolution under multi-step thermo-mechanical schedules,
2. forward field reconstruction under strong spatial heterogeneity, 
3. inverse identification under classic isothermal relaxation,
4. inverse identification under non-isothermal schedules to calibrate time–temperature coupling.

---

