function kin = planar_kinematics(Q, p)
%PLANAR_KINEMATICS Rigid-body COM positions and velocity Jacobians.
% Input Q = [xb; yb; psi; q1; q2] [m,m,rad,rad,rad].
% Positions and translational Jacobians are expressed in the inertial frame.

assert(numel(Q) == 5, 'FFSM:Dimension', 'Q must contain five coordinates.');
Q = Q(:);
xb = Q(1); yb = Q(2); psi = Q(3); q1 = Q(4); q2 = Q(5);
S = [0 -1; 1 0];
e = [1; 0];
R0 = rot2(psi);
R1 = rot2(psi + q1);
R2 = rot2(psi + q1 + q2);
r = [xb; yb];

d0 = R0 * (p.a0 * e);
d1c = R1 * (p.c1 * e);
d1 = R1 * (p.l1 * e);
d2c = R2 * (p.c2 * e);
d2 = R2 * (p.l2 * e);

kin.positions = [r, r+d0+d1c, r+d0+d1+d2c];
kin.Jv = cell(1,3);
kin.Jw = cell(1,3);
kin.Jv{1} = [eye(2), zeros(2,3)];
kin.Jv{2} = [eye(2), S*d0+S*d1c, S*d1c, zeros(2,1)];
kin.Jv{3} = [eye(2), S*d0+S*d1+S*d2c, S*d1+S*d2c, S*d2c];
kin.Jw{1} = [0 0 1 0 0];
kin.Jw{2} = [0 0 1 1 0];
kin.Jw{3} = [0 0 1 1 1];
kin.masses = [p.m0 p.m1 p.m2];
kin.inertias = [p.J0 p.J1 p.J2];
kin.names = {'base','link1','link2'};
kin.shoulder = r + d0;
kin.elbow = r + d0 + d1;
kin.end_effector = r + d0 + d1 + d2;

if p.payload.enabled
    kin.positions(:,end+1) = kin.end_effector;
    kin.Jv{end+1} = [eye(2), S*d0+S*d1+S*d2, S*d1+S*d2, S*d2];
    kin.Jw{end+1} = [0 0 1 1 1];
    kin.masses(end+1) = p.payload.mass;
    kin.inertias(end+1) = p.payload.J;
    kin.names{end+1} = 'payload';
end
end

function R = rot2(angle)
R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
end
