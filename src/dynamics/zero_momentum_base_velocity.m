function bdot = zero_momentum_base_velocity(Q, qdot, p)
%ZERO_MOMENTUM_BASE_VELOCITY Enforce P=[0;0], L_O=0.
% This relation is valid only for the unforced zero-momentum level.

H = full_mass_matrix(Q, p);
A = H(1:3,1:3);
Bc = H(1:3,4:5);
bdot = -(A \ (Bc * qdot(:)));
end
