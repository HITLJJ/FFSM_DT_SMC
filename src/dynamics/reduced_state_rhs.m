function xdot = reduced_state_rhs(t, x, p, torque_fun)
%REDUCED_STATE_RHS Seven-state zero-momentum shape dynamics and reconstruction.
% x=[b(3);q(2);qdot(2)]. No independent base input is accepted.

b = x(1:3); q = x(4:5); qdot = x(6:7);
Q = [b; q];
bdot = zero_momentum_base_velocity(Q, qdot, p);
fullState = [Q; bdot; qdot];
tau = torque_fun(t, fullState);
Mr = reduced_mass_matrix(q, p);
C = reduced_coriolis_matrix(q, qdot, p);
qddot = Mr \ (tau(:) - C*qdot);
xdot = [bdot; qdot; qddot];
end
