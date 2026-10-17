# FFSM_DT_SMC

This repository develops a planar 3+2 free-floating robot from validated
dynamics through certified uncertainty bounds to an exploratory
Plant-Digital-Twin-SMC closed loop.

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
run_stage1_tests;                 % validated full/reduced dynamics
run_stage2_tests;                 % certified reference-torque baseline
run_stage3a_tests;                % certified Taylor-Hessian bounds
run_closed_loop_tests;            % Plant-DT-SMC interface/regression tests
summary = run_closed_loop_comparison; % MAT, CSV, FIG and PNG
```

Closed-loop output is written to `results/closed_loop/`; earlier stage outputs
remain in their respective result directories. The optional Stage-1 independent
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
- `src/uncertainty/`: certified and exploratory bound interfaces
- `src/controllers/`: fixed, adaptive and bound-scheduled SMC
- `src/digital_twin/`: isolated synchronization and exploratory online updates

Stage 2 and Stage 3A certificates remain distinct from exploratory closed-loop
methods. The local Taylor candidate and online parameter-set update are never
reported as certified. See `docs/closed_loop_research_report.md` for results,
limitations and the strict/exploratory boundary.
