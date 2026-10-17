function blocks=affine_mass_block_derivatives(q,theta,geometry)
%AFFINE_MASS_BLOCK_DERIVATIVES Parameter coefficients of H at fixed q.
% q is 2x1 [rad]. theta=[m0,m1,m2,J0,J1,J2] uses kg and kg*m^2.
% The full inertia is affine in theta, so the returned A_j, B_j and D_j
% are analytic parameter derivatives. This floating implementation is for
% derivative cross-validation, not for interval certification.

theta=theta(:); q=q(:);
p=inertial_parameters_to_model(theta,geometry);
Q=[zeros(3,1,'like',q);q];
H=full_mass_matrix(Q,p);
blocks.A=H(1:3,1:3); blocks.B=H(1:3,4:5); blocks.D=H(4:5,4:5);
blocks.Aj=zeros(3,3,6,'like',H);
blocks.Bj=zeros(3,2,6,'like',H);
blocks.Dj=zeros(2,2,6,'like',H);
for j=1:6
    shifted=theta; shifted(j)=shifted(j)+1;
    Hj=full_mass_matrix(Q,inertial_parameters_to_model(shifted,geometry))-H;
    blocks.Aj(:,:,j)=Hj(1:3,1:3);
    blocks.Bj(:,:,j)=Hj(1:3,4:5);
    blocks.Dj(:,:,j)=Hj(4:5,4:5);
end
end
