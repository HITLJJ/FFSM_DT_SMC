function result = simulate_stage1_case(caseName, cfg)
%SIMULATE_STAGE1_CASE Run full and zero-momentum reduced dynamics.

p = cfg.parameters;
b0 = [0; 0; 0.10];
q0 = [0.40; -0.70];
switch upper(caseName)
    case 'A'
        qdot0 = [0;0];
    case 'B'
        qdot0 = [0.25;-0.18];
    case 'C'
        qdot0 = [0;0];
    otherwise
        error('FFSM:UnknownCase', 'Unknown experiment case %s.', caseName);
end
Q0 = [b0;q0];
bdot0 = zero_momentum_base_velocity(Q0, qdot0, p);
x0 = [Q0;bdot0;qdot0];
xr0 = [b0;q0;qdot0];
tspan = (0:cfg.output_dt:cfg.t_final).';
torqueFun = @(t,x) stage1_torque(caseName,t,x);
rhs = @(t,x) full_state_rhs(t,x,p,torqueFun);
% Integrate actuator work with the same adaptive solver rather than applying
% trapezoidal quadrature only at requested plot/output times. This keeps the
% work-energy acceptance metric independent of output sampling density.
rhsWithWork = @(t,z) [rhs(t,z(1:10)); ...
    z(9:10).'*torqueFun(t,z(1:10))];
rhsReduced = @(t,x) reduced_state_rhs(t,x,p,torqueFun);
solver = str2func(cfg.solver);

strictOptions = odeset('RelTol',cfg.strict.RelTol, ...
    'AbsTol',cfg.strict.AbsTol,'MaxStep',cfg.strict.MaxStep);
coarseOptions = odeset('RelTol',cfg.coarse.RelTol, ...
    'AbsTol',cfg.coarse.AbsTol,'MaxStep',cfg.coarse.MaxStep);

[tStrict,zStrict] = solver(rhsWithWork,tspan,[x0;0],strictOptions);
[tCoarse,zCoarse] = solver(rhsWithWork,tspan,[x0;0],coarseOptions);
[tReduced,xReduced] = solver(rhsReduced,tspan,xr0,strictOptions);
xStrict = zStrict(:,1:10);
xCoarse = zCoarse(:,1:10);

result.case_name = upper(caseName);
result.t = tStrict;
result.x = xStrict;
result.actuator_work = zStrict(:,11);
result.t_reduced = tReduced;
result.x_reduced = xReduced;
result.convergence.max_state_difference = max(vecnorm((xStrict-xCoarse).',2,1));
result.convergence.final_state_difference = norm(xStrict(end,:)-xCoarse(end,:));
result.convergence.final_work_difference = abs(zStrict(end,11)-zCoarse(end,11));
result.solver.strict = cfg.strict;
result.solver.coarse = cfg.coarse;
result.diagnostics = compute_stage1_diagnostics(result,p,torqueFun);
end
