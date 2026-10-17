# FFSM_DT_SMC

Stage 1 implements and validates the planar 3+2 free-floating robot dynamics

```text
Q = [x_b, y_b, psi, q_1, q_2]^T,
B = [zeros(3,2); eye(2)].
```

Only the two relative joints are actuated. There is no base thrust or base-yaw
actuator. The baseline has three rigid bodies (base and two links), no gravity,
friction, contact or payload, and uses SI units.

## MATLAB R2025b quick start

From the repository root:

```matlab
addpath(genpath(pwd));
results = run_stage1_tests;       % automated acceptance tests
summary = run_stage1_validation;  % experiments A-C, MAT, FIG, PNG and AVI
```

The experiment output is written to `results/stage1/`. Numerical settings and
acceptance thresholds are in `config/stage1_config.m`. The optional independent
symbolic cross-check is:

```matlab
report = symbolic_consistency_check;
```

That check requires Symbolic Math Toolbox. The runtime dynamics do not require
symbolic operations and never use `inv()`.

## Layout

- `direct_dynamics-1.0.0/`: untouched public reference implementation
- `manuscript/`: theory source and PDF
- `config/`: parameters, solver settings and thresholds
- `src/dynamics/`: complete and zero-momentum reduced dynamics
- `tests/`: automated dynamics acceptance tests
- `experiments/`: experiments A-C, diagnostics, plots and animation
- `results/`: generated reproducibility artifacts
- `docs/`: stage reports
- `references/`: additional reference material
- `src/uncertainty/`, `src/controllers/`, `src/digital_twin/`: reserved for later stages

Stage 1 deliberately does not implement sliding-mode control, certified
uncertainty bounds, online identification or digital-twin model updates.
