function blocks = partition_dynamics(Q, Qdot, p)
%PARTITION_DYNAMICS Passive-base/active-joint partition of H and h.
H = full_mass_matrix(Q, p);
h = full_nonlinear_vector(Q, Qdot, p);
blocks.H = H;
blocks.h = h;
blocks.A = H(1:3,1:3);
blocks.Bc = H(1:3,4:5);
blocks.D = H(4:5,4:5);
blocks.hb = h(1:3);
blocks.hq = h(4:5);
blocks.Mr = blocks.D - blocks.Bc.' * (blocks.A \ blocks.Bc);
blocks.hr = blocks.hq - blocks.Bc.' * (blocks.A \ blocks.hb);
end
