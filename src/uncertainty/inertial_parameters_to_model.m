function p = inertial_parameters_to_model(theta, geometry)
%INERTIAL_PARAMETERS_TO_MODEL Apply inertial parameters to fixed geometry.
% theta=[m0;m1;m2;J0;J1;J2], with kg and kg*m^2 units. geometry follows
% default_parameters and is not mutated. The Stage-2 baseline has no payload.

theta = theta(:);
assert(numel(theta)==6,'FFSM:Dimension','theta must have six elements.');
assert(all(isfinite(real(theta))) && all(real(theta)>0), ...
    'FFSM:InvalidInertialParameters','All masses and inertias must be positive finite values.');
if ~isreal(theta)
    assert(all(abs(imag(theta)) <= 1e-6*max(1,abs(real(theta)))), ...
        'FFSM:InvalidInertialParameters','Complex perturbations are too large.');
end
p = geometry;
p.m0=theta(1); p.m1=theta(2); p.m2=theta(3);
p.J0=theta(4); p.J1=theta(5); p.J2=theta(6);
p.payload.enabled=false; p.payload.mass=0; p.payload.J=0;
end
