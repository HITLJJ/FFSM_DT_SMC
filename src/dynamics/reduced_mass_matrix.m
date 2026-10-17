function Mr = reduced_mass_matrix(q, p)
%REDUCED_MASS_MATRIX Zero-momentum two-joint Schur-complement inertia.
% The SE(2)-invariant shape metric depends only on q=[q1;q2].

q = q(:);
assert(numel(q) == 2, 'FFSM:Dimension', 'q must contain two joint angles.');
Q = [zeros(3,1,'like',q); q];
H = full_mass_matrix(Q, p);
A = H(1:3,1:3);
Bc = H(1:3,4:5);
D = H(4:5,4:5);
Mr = D - Bc.' * (A \ Bc);
end
