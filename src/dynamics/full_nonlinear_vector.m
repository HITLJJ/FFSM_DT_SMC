function h = full_nonlinear_vector(Q, Qdot, p)
%FULL_NONLINEAR_VECTOR Exact velocity-quadratic Euler-Lagrange term.
% h is computed from independent rigid-body acceleration-bias expressions,
% not by differentiating the implemented mass-matrix entries.

assert(numel(Q) == 5 && numel(Qdot) == 5, 'FFSM:Dimension', ...
    'Q and Qdot must each contain five elements.');
Q = Q(:); Qdot = Qdot(:);
psi = Q(3); q1 = Q(4); q2 = Q(5);
dpsi = Qdot(3); dq1 = Qdot(4); dq2 = Qdot(5);
e = [1;0];
R0 = rot2(psi); R1 = rot2(psi+q1); R2 = rot2(psi+q1+q2);

a0bias = zeros(2,1, 'like', Qdot);
a1bias = -R0*(p.a0*e)*dpsi^2 ...
    -R1*(p.c1*e)*(dpsi+dq1)^2;
a2bias = -R0*(p.a0*e)*dpsi^2 ...
    -R1*(p.l1*e)*(dpsi+dq1)^2 ...
    -R2*(p.c2*e)*(dpsi+dq1+dq2)^2;

kin = planar_kinematics(Q, p);
biases = {a0bias, a1bias, a2bias};
if p.payload.enabled
    biases{end+1} = -R0*(p.a0*e)*dpsi^2 ...
        -R1*(p.l1*e)*(dpsi+dq1)^2 ...
        -R2*(p.l2*e)*(dpsi+dq1+dq2)^2;
end

h = zeros(5,1, 'like', Qdot);
for i = 1:numel(kin.masses)
    h = h + kin.masses(i) * kin.Jv{i}.' * biases{i};
end
end

function R = rot2(angle)
R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
end
