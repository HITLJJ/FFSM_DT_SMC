# Stage 1: dynamics implementation and validation report

## Status and scope

This stage implements the ideal planar 3+2 free-floating rigid-body system with
five generalized coordinates and two internal actuator inputs. Sliding-mode
control, uncertainty certification, online parameter identification and digital
twin updates are intentionally excluded.

The implementation was accepted on **2026-10-17** using MATLAB
`25.2.0.2998904 (R2025b)` on Windows 11. MATLAB and Symbolic Math Toolbox
licenses were active. The installed relevant toolboxes include Symbolic Math,
Control System, Robust Control, Robotics System, Optimization, Global
Optimization, Parallel Computing, Simscape and Simscape Multibody.

All 10 automated tests passed after the validation-instrument corrections
documented below. The independent symbolic check also passed. The numerical
results support proceeding to Stage 2 within the same ideal planar 3+2 model
scope; they are regression evidence, not a replacement for the analytical
proofs in the theory document.

## Source material inspected

- Theory norm: `manuscript/FFSM_3plus2_Theory_CN.pdf` (13 pages). The requested
  name contained `(4)`, while the file actually present does not. Its content
  matches the requested 3+2 theory: COM kinematics, Jacobian inertia, full
  Euler-Lagrange dynamics, momentum map, Schur complement and zero-momentum
  shape reduction.
- Public reference: `direct_dynamics-1.0.0/Openloop_FFSM_2Links.m` by S. R.
  Nekoo. The original file and its citation/license notice are unchanged.

## Coordinates, geometry and parameters

The inertial-frame generalized coordinate is

\[
Q=[x_b,y_b,\psi,q_1,q_2]^T.
\]

`(x_b,y_b)` is the base center of mass, `psi` is base yaw, and `q1,q2` are
relative joint angles. The shoulder lies `a0` metres along the base-frame
positive x-axis. The model uses the theory document's positions

\[
p_0=r,\quad p_1=r+R(\psi)a_0e+R(\psi+q_1)c_1e,
\]

\[
p_2=r+R(\psi)a_0e+R(\psi+q_1)\ell_1e
       +R(\psi+q_1+q_2)c_2e.
\]

Baseline values are `m0=5 kg`, `m1=m2=2 kg`, `l1=l2=2 m`, `c1=c2=1 m`,
`a0=0.9 m`, with uniform-link inertias and the rectangular-base inertia listed
in `src/dynamics/default_parameters.m`. Geometry is treated as known. The
uncertain inertial parameter vector reserved for later stages is
`theta=[m0,m1,m2,J0,J1,J2]^T`.

The payload configuration exists but is disabled (`mass=0`) in the baseline.
When enabled it is added exactly once at the link-2 tip as either a point mass
or a planar rigid body.

## Inputs and degrees of freedom

The system has five generalized degrees of freedom and only two independent
inputs:

\[
B\tau=\begin{bmatrix}0_{3\times2}\\I_2\end{bmatrix}\tau.
\]

No API in Stage 1 accepts an independent base force or base-yaw torque.

## Changes relative to the public code

1. The public code's two-component shoulder offset `(ex,ey)` is replaced by the
   theory document's scalar `a0` along the base positive x-axis. This is an
   explicit coordinate-model choice, not an elementwise formula patch.
2. The public default payload `mp=1 kg` is disabled; the baseline is exactly
   three rigid bodies. Optional payload bookkeeping prevents duplicate mass.
3. The handwritten expanded inertia and nonlinear terms are replaced by COM
   Jacobian assembly and independent rigid-body acceleration-bias assembly.
4. Explicit `inv()` is eliminated in favor of matrix left division.
5. Explicit Euler integration is replaced by `ode113` with strict and coarse
   tolerance sets. Both runs share requested output times for convergence
   comparison.
6. Nonzero initial joint rates are paired with base rates computed from the
   zero-momentum connection; the public code's unrelated base rates are not
   reused.
7. Momentum, angular momentum, direct kinetic energy, work, model consistency,
   minimum eigenvalue and reconstruction diagnostics are recorded.

## Dynamics implementation

`full_mass_matrix` assembles

\[
H=\sum_i(m_iJ_{v,i}^TJ_{v,i}+J_i j_{\omega,i}^Tj_{\omega,i}).
\]

`full_nonlinear_vector` independently forms the generalized inertial bias from
`sum(Jv_i' * m_i * a_bias_i)`. This is exactly the velocity-quadratic
Euler-Lagrange term for the stated rigid-body kinematics; the angular Jacobians
are constant, so their acceleration bias is zero.

The passive/active partition is

\[
H=\begin{bmatrix}A&B_c\\B_c^T&D\end{bmatrix},\qquad
M_r=D-B_c^TA^{-1}B_c.
\]

All code evaluates products with `A\Bc` or linear solves. The zero-momentum
connection is `bdot=-(A\(Bc*qdot))`. The reduced Coriolis matrix is constructed
from the Christoffel relation for `Mr`; complex-step derivatives give stable
numeric derivatives without symbolic work inside ODE calls.

Two independent validation paths are supplied:

- direct rigid-body kinetic energy and acceleration-bias dynamics versus the
  mass-matrix representation;
- an optional Symbolic Math Toolbox script that derives `H` and `h` from the
  kinetic energy, plus a central finite-difference Christoffel test.

## Experiments

- **A - completely stationary:** zero generalized velocity and zero torque.
- **B - zero-momentum free motion:** `qdot(0)=[0.25,-0.18] rad/s`; base velocity
  is reconstructed from the zero-momentum connection; torque is zero.
- **C - internal torque driven:** zero initial velocity and deterministic smooth
  bounded joint torques. Momentum must remain constant, while kinetic-energy
  change must match integrated actuator power.

The default solver is `ode113`, duration `10 s`, output interval `0.01 s`.
Strict settings are `RelTol=1e-10`, `AbsTol=1e-12`, `MaxStep=0.01 s`; coarse
settings are `RelTol=1e-7`, `AbsTol=1e-9`, `MaxStep=0.05 s`.

Actuator work is integrated as an additional scalar ODE state
`Wdot=qdot'*tau` using the same adaptive solver. Sampled power is still retained
for plotting. This prevents the work-energy acceptance result from depending on
the plot/output sampling interval.

## Acceptance table

Thresholds are absolute engineering regression tolerances for the baseline SI
parameter scale. They are not mathematical proofs. They should be revisited if
the system is rescaled or parameters change by orders of magnitude.

| Test | Metric | Threshold | MATLAB R2025b result | Status |
|---|---:|---:|---:|---|
| Inertia symmetry | max `norm(H-H',2)` | `1e-12` | `0` | PASS |
| Full inertia SPD | sampled min `eig(H)` | `>1e-8` | `0.347312883969385` | PASS |
| Reduced inertia SPD | sampled min `eig(Mr)` | `>1e-8` | `0.3651944963716924` | PASS |
| Full/reduced acceleration | max `norm(ddq_full-ddq_red)` | `<1e-8 rad/s^2` | `1.233896894832203e-14 rad/s^2` | PASS |
| Linear momentum | max `norm(P-P0)` | `<1e-7 kg m/s` | `9.461842752807287e-15 kg m/s` | PASS |
| Angular momentum | max `abs(LO-LO0)` | `<1e-7 kg m^2/s` | `3.145400606641147e-14 kg m^2/s` | PASS |
| Free-motion energy | max `abs(T-T0)` | `<1e-7 J` | `2.567390744445674e-16 J` | PASS |
| Driven work-energy | max `abs((T-T0)-W)` | `<1e-6 J` | `1.398881011027697e-14 J` | PASS |
| Nonlinear finite difference | `norm(h-h_fd)` | `<1e-6` generalized-force units | `3.501615323840983e-9` | PASS |
| Zero-momentum reconstruction | max `norm(bdot-bdot_reconstructed)` | `<1e-8` mixed SI rates | `1.597763194179639e-15` | PASS |
| Solver convergence | max strict/coarse state 2-norm | `<1e-5` mixed SI state units | `4.32652216294636e-9` | PASS |

The symmetry threshold is near double-precision roundoff because `H` is a sum
of symmetric Gram matrices. Conservation and work-energy thresholds allow for
adaptive integration error. The finite-difference threshold is looser because
it includes truncation and cancellation error.

Additional independent checks passed:

- symbolic-versus-numeric mass-matrix error: `2.0940e-15`;
- symbolic-versus-numeric nonlinear-vector error: `7.8505e-17`;
- direct rigid-body versus matrix kinetic-energy error: `2.131628207280301e-14 J`;
- optional payload single-counting error: `3.951940768217079e-15`;
- experiment-B initial zero-momentum residuals: `norm(P)=8.122e-17 kg m/s`,
  `abs(LO)=5.065e-16 kg m^2/s`.

### Per-experiment extrema

| Metric | A: stationary | B: free motion | C: torque driven |
|---|---:|---:|---:|
| max linear-momentum drift [kg m/s] | `0` | `9.461842752807287e-15` | `5.092297236337201e-15` |
| max angular-momentum drift [kg m^2/s] | `0` | `3.145400606641147e-14` | `1.898481372109018e-14` |
| max kinetic-energy change [J] | `0` | `2.567390744445674e-16` | `0.2279177964556557` |
| max work-energy residual [J] | `0` | `2.567390744445674e-16` | `1.398881011027697e-14` |
| max full/reduced acceleration error [rad/s^2] | `0` | `2.120907854283028e-16` | `3.436427221258261e-15` |
| min eigenvalue of `H` | `0.5485273487169101` | `0.5485273487169101` | `0.4529073994680357` |
| min eigenvalue of `Mr` | `0.5960406296276324` | `0.5960406296276324` | `0.4677631020943784` |
| strict/coarse max state difference | `0` | `6.336574504657514e-12` | `4.32652216294636e-9` |

Experiment C is not expected to conserve kinetic energy: the observed
`0.2279177964556557 J` change is accounted for by internal actuator work to
within `1.40e-14 J`.

## Validation anomalies and resolutions

The initial acceptance run is retained here for transparency.

1. The first experiment-C unit test reported a `2.638259698041837e-6 J`
   work-energy residual against the `1e-6 J` threshold. A 10-second output-step
   study gave `2.8638e-5`, `7.1588e-6`, `1.7896e-6`, `4.4742e-7` and
   `1.1185e-7 J` for output intervals `0.04`, `0.02`, `0.01`, `0.005` and
   `0.0025 s`, respectively. The second-order decrease identified sampled
   trapezoidal quadrature as the cause. Work was therefore added as an ODE
   state. Neither dynamics nor the `1e-6 J` threshold was changed.
2. The first symbolic check failed before evaluation because R2025b reserves
   `psi` as a function name inside `syms`. The test variable was renamed to
   `psi_angle`; the resulting independent symbolic errors were then at machine
   precision. No dynamics expression was changed.
3. MATLAB startup emits a Miniconda `requests` dependency warning. It did not
   affect MATLAB exit status, numerical tests, figures or video generation, but
   the Python environment should be cleaned separately before relying on
   MATLAB-Python integration.

## Outputs and reproducibility

`run_stage1_validation` saves, for each experiment:

- raw trajectory and diagnostics in MAT format;
- native FIG and PNG state, momentum, energy and consistency plots;
- experiment C planar animation as Motion JPEG AVI;
- a machine-readable summary MAT file and a plain-text numerical summary.

All requested artifacts were generated successfully in `results/stage1/`.
The experiment-C AVI was reopened using `VideoReader`: duration
`11.13333333333333 s`, frame rate `30 fps`, resolution `1394x840`. Its playback
duration differs from the simulated 10 seconds because trajectory frames were
subsampled and encoded at a fixed presentation frame rate.

## Unresolved items

1. If the physical platform uses a nonzero lateral shoulder offset, extend the
   theory and implementation consistently instead of silently restoring `ey`.
2. Validate real hardware mass properties and reachable geometry before using
   the baseline numerical values in an experiment.
3. Stage 1 verifies the nominal ideal model only. It does not certify parameter
   uncertainty sets, unmodelled forces, sensor noise or physical hardware.
4. The random-state and trajectory tests are numerical evidence only; global
   symmetry/positive-definiteness conclusions remain those of the analytical
   model, not Monte Carlo proofs.

## Stage 2 decision

**Recommended, with scope restrictions.** Stage 1 acceptance is complete for
the nominal ideal planar 3+2 free-floating system, so work may proceed to Stage
2: construction and validation of strict uncertainty bounds. Stage 2 must keep
the verified plant model separate from the controller's nominal model, must not
treat sampling as a proof of a strict bound, and must not introduce independent
base actuation. This recommendation does not authorize Stage 3 sliding-mode
control or later digital-twin updates.
