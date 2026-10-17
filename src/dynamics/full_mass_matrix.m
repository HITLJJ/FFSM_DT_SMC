function H = full_mass_matrix(Q, p)
%FULL_MASS_MATRIX Complete 5-by-5 generalized inertia matrix [mixed SI units].
% Constructed independently from rigid-body COM velocity Jacobians.

kin = planar_kinematics(Q, p);
H = zeros(5,5, 'like', Q);
for i = 1:numel(kin.masses)
    Jv = kin.Jv{i};
    Jw = kin.Jw{i};
    H = H + kin.masses(i) * (Jv.' * Jv) + ...
        kin.inertias(i) * (Jw.' * Jw);
end
% Do not force symmetry numerically: tests must be able to detect asymmetry.
end
