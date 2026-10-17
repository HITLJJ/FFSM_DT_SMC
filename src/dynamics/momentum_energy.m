function out = momentum_energy(Q, Qdot, p)
%MOMENTUM_ENERGY Independent total momentum and kinetic energy calculation.
% P is inertial linear momentum [kg*m/s]; LO is angular momentum about the
% fixed inertial origin [kg*m^2/s]; T is kinetic energy [J].

Q = Q(:); Qdot = Qdot(:);
kin = planar_kinematics(Q, p);
P = zeros(2,1);
LO = 0;
T = 0;
bodyVelocity = cell(1,numel(kin.masses));
for i = 1:numel(kin.masses)
    vi = kin.Jv{i} * Qdot;
    wi = kin.Jw{i} * Qdot;
    pi = kin.positions(:,i);
    P = P + kin.masses(i)*vi;
    LO = LO + pi(1)*(kin.masses(i)*vi(2)) ...
        - pi(2)*(kin.masses(i)*vi(1)) + kin.inertias(i)*wi;
    T = T + 0.5*kin.masses(i)*(vi.'*vi) + 0.5*kin.inertias(i)*wi^2;
    bodyVelocity{i} = vi;
end
out.P = P;
out.LO = LO;
out.T = T;
out.body_velocity = bodyVelocity;
end
