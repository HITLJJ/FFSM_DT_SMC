function [F,details] = reference_torque_map(z, theta, geometry)
%REFERENCE_TORQUE_MAP Control-related zero-momentum reference torque [N*m].
% z is a struct with 2x1 fields:
%   q [rad], qdot [rad/s], qr_dot [rad/s], qr_ddot [rad/s^2].
% theta=[m0;m1;m2;J0;J1;J2] uses kg and kg*m^2. Geometry is fixed in SI.
% Output F=Mr*qr_ddot+Cr*qr_dot is 2x1. It never uses a simulated "true"
% parameter and must not be replaced by Mr*qr_ddot+Cr*qdot.

required={'q','qdot','qr_dot','qr_ddot'};
assert(isstruct(z) && all(isfield(z,required)),'FFSM:InvalidReferenceState', ...
    'z must contain q, qdot, qr_dot and qr_ddot.');
for k=1:numel(required)
    value=z.(required{k});
    assert(isnumeric(value) && numel(value)==2 && all(isfinite(value)), ...
        'FFSM:InvalidReferenceState','Every z field must be a finite 2x1 vector.');
end
p=inertial_parameters_to_model(theta,geometry);
q=z.q(:); qdot=z.qdot(:); qrDot=z.qr_dot(:); qrDDot=z.qr_ddot(:);
Mr=reduced_mass_matrix(q,p);
C=reduced_coriolis_matrix(q,qdot,p);
F=Mr*qrDDot+C*qrDot;
if nargout>1
    details.Mr=Mr; details.C=C;
    details.inertia_component=Mr*qrDDot;
    details.coriolis_component=C*qrDot;
end
end
