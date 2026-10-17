function Qddot = full_acceleration(Q, Qdot, tau, p)
%FULL_ACCELERATION Solve H*Qddot+h=B*tau without explicit inversion.
H = full_mass_matrix(Q, p);
h = full_nonlinear_vector(Q, Qdot, p);
Qddot = H \ (input_matrix()*tau(:) - h);
end
