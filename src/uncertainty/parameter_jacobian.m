function J = parameter_jacobian(z, theta, geometry)
%PARAMETER_JACOBIAN Central-difference dF/dtheta reference calculation.
% Output is 2x6 with mixed torque-per-parameter units. This is a numerical
% sensitivity diagnostic, not a certified interval derivative.

theta=theta(:); J=zeros(2,6);
for k=1:6
    h=2e-5*max(abs(theta(k)),1);
    d=zeros(6,1); d(k)=h;
    J(:,k)=(reference_torque_map(z,theta+d,geometry)- ...
        reference_torque_map(z,theta-d,geometry))/(2*h);
end
end
