function report = symbolic_consistency_check()
%SYMBOLIC_CONSISTENCY_CHECK Optional independent Symbolic Math Toolbox check.
% This is intentionally separate from the simulation implementation. It
% constructs H and h from a symbolic kinetic energy, then compares them at
% a deterministic state with the numerical rigid-body implementation.

assert(license('test','Symbolic_Toolbox') && ~isempty(ver('symbolic')), ...
    'FFSM:SymbolicUnavailable','Symbolic Math Toolbox is required.');
p = default_parameters();
xb = sym('xb','real'); yb = sym('yb','real');
psiSym = sym('psi_angle','real');
q1 = sym('q1','real'); q2 = sym('q2','real');
dxb = sym('dxb','real'); dyb = sym('dyb','real');
dpsi = sym('dpsi','real'); dq1 = sym('dq1','real'); dq2 = sym('dq2','real');
Q = [xb;yb;psiSym;q1;q2]; dQ = [dxb;dyb;dpsi;dq1;dq2];
R = @(a) [cos(a) -sin(a);sin(a) cos(a)]; e=[1;0];
positions = {Q(1:2), Q(1:2)+R(psiSym)*(p.a0*e)+R(psiSym+q1)*(p.c1*e), ...
    Q(1:2)+R(psiSym)*(p.a0*e)+R(psiSym+q1)*(p.l1*e)+R(psiSym+q1+q2)*(p.c2*e)};
angles = {psiSym,psiSym+q1,psiSym+q1+q2};
masses = [p.m0 p.m1 p.m2]; inertias = [p.J0 p.J1 p.J2];
T = sym(0);
for i = 1:3
    v = jacobian(positions{i},Q)*dQ;
    w = jacobian(angles{i},Q)*dQ;
    T = T + masses(i)*(v.'*v)/2 + inertias(i)*w^2/2;
end
Hsym = simplify(hessian(T,dQ));
dTdQdot = jacobian(T,dQ).';
hSym = simplify(jacobian(dTdQdot,Q)*dQ-jacobian(T,Q).');

Qn = [0.2;-0.4;0.3;0.5;-0.7]; dQn = [-0.1;0.2;-0.05;0.3;-0.25];
vars = [Q;dQ]; values = [Qn;dQn];
Hcheck = double(subs(Hsym,vars,values));
hCheck = double(subs(hSym,vars,values));
report.mass_error = norm(Hcheck-full_mass_matrix(Qn,p),2);
report.nonlinear_error = norm(hCheck-full_nonlinear_vector(Qn,dQn,p),2);
report.pass = report.mass_error < 1e-10 && report.nonlinear_error < 1e-10;
assert(report.pass,'FFSM:SymbolicMismatch','Symbolic and numerical dynamics differ.');
end
