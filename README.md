# Reduced-Time Recursive PINN for Thermo-Viscoelastic SMPCs

**Author:** Haobo Wang, Harbin Institute of Technology (HIT), Harbin, China  
**Email:** wanghaobohit@163.com

## Overview

This repository supports a mechanics-oriented study on thermo-viscoelastic shape memory polymer composites (SMPCs). The central problem is forward field reconstruction and inverse constitutive identification under history-dependent thermo-mechanical loading.

The proposed framework, referred to as a reduced-time recursive physics-informed neural network (RT-RPINN), embeds the reduced-time Maxwell state update into the learning process. The neural network represents the displacement field, while strain, stress, internal variables, equilibrium residuals, boundary tractions, and global reaction quantities are generated through constitutive-consistent recursive updates.

![Overview](Overview.png)

## Mechanical Problem

SMPCs exhibit coupled responses controlled by:

- anisotropic stiffness and thermal expansion,
- viscoelastic relaxation and internal material memory,
- temperature-dependent WLF/Arrhenius shift behavior,
- constrained-stage reaction forces and recovery-stage deformation.

The main inverse task is to identify thermo-viscoelastic constitutive parameters from limited observations while preserving consistency with the governing thermo-mechanical behavior. Particular attention is paid to parameter observability, loading-path design, and robustness across random seeds.

## Method

The framework is used for two related tasks:

- **Forward reconstruction:** recover spatiotemporal displacement fields and derived strain/stress responses from sparse observations and fixed material parameters.
- **Inverse identification:** recover selected constitutive parameters from sparse displacement observations, optional strain information, and optional reaction-force histories.

Two temporal representations are compared:

- **RT-RPINN:** a standard spatiotemporal neural backbone coupled to a reduced-time recursive constitutive module.
- **LSTM-based RT-RPINN:** a sequence-aware temporal backbone used to test whether neural temporal memory improves forward reconstruction or inverse identification.

Reaction-force supervision is treated as an additional observable, not as a universally beneficial constraint. Its usefulness depends on whether it provides independent stage-specific information beyond displacement and strain observations.

## Benchmark Examples

- **EX1:** Forward reconstruction of a non-isothermal shape-memory cycle.
- **EX2:** Inverse identification under isothermal ramp-hold stress-relaxation loading.
- **EX3:** Forward reconstruction in a heterogeneous stress/strain field case.
- **EX4:** Inverse identification under non-isothermal loading, including multi-seed RF/no-RF evaluation.

Together, these examples examine displacement-field accuracy, derived strain/stress sensitivity, mechanical relaxation identification, thermal-shift identifiability, and the role of reaction-force observables.

## Repository Structure

```text
EX1/                 Forward non-isothermal shape-memory benchmark
EX2/                 Isothermal inverse stress-relaxation benchmark
EX3/                 Heterogeneous forward reconstruction benchmark
EX4/                 Non-isothermal inverse benchmark and multi-seed analysis
requirements.txt     Python dependency list
datasavebystep.py    Utility for step-wise data export and organization
umat.for             Fortran UMAT-side constitutive implementation
materialtable.txt    Material parameters used by the constitutive setup
Overview.png         Schematic overview of the framework
```

## Notes

The examples are designed as controlled computational benchmarks. Reference responses are generated from finite-element/constitutive simulations so that reconstruction errors, parameter bias, and identifiability trends can be evaluated against known material parameters.

When using or extending the code, please report:

- the observation density and time-frame sampling strategy,
- whether reaction-force supervision is used,
- the random seeds used for inverse identification,
- the trainable constitutive parameter subset,
- the material parameters and temperature-loading programme.
