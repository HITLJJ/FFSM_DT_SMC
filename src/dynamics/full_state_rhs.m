function xdot = full_state_rhs(t, x, p, torque_fun)
%FULL_STATE_RHS Complete ten-state free-floating dynamics.
% x=[Q;Qdot], torque_fun(t,x) returns two internal joint torques [N*m].

Q = x(1:5); Qdot = x(6:10);
tau = torque_fun(t, x);
H = full_mass_matrix(Q, p);
h = full_nonlinear_vector(Q, Qdot, p);
Qddot = H \ (input_matrix()*tau(:) - h);
xdot = [Qdot; Qddot];
end
