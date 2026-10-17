function qddot = reduced_acceleration(q, qdot, tau, p)
%REDUCED_ACCELERATION Solve zero-momentum shape dynamics.
Mr = reduced_mass_matrix(q, p);
C = reduced_coriolis_matrix(q, qdot, p);
qddot = Mr \ (tau(:) - C*qdot(:));
end
