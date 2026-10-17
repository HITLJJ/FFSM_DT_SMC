function p = default_parameters()
%DEFAULT_PARAMETERS Baseline planar 3+2 free-floating robot parameters.
% Coordinates follow FFSM_3plus2_Theory_CN.pdf. The base COM is Q(1:2),
% the shoulder is a0 metres along the base +x axis, and q1/q2 are relative
% revolute angles. All quantities use SI units.

p.m0 = 5.0;                  % base mass [kg]
p.m1 = 2.0;                  % link 1 mass [kg]
p.m2 = 2.0;                  % link 2 mass [kg]
p.l1 = 2.0;                  % joint 1 to joint 2 [m]
p.l2 = 2.0;                  % joint 2 to end effector [m]
p.c1 = p.l1 / 2;             % joint 1 to link-1 COM [m]
p.c2 = p.l2 / 2;             % joint 2 to link-2 COM [m]
p.a0 = 0.9;                  % base COM to shoulder along base +x [m]
p.base_half_length = 1.0;    % animation geometry [m]
p.base_half_width = 0.5;     % animation geometry [m]
p.J0 = p.m0 * ((2*p.base_half_length)^2 + ...
    (2*p.base_half_width)^2) / 12;                  % [kg*m^2]
p.J1 = p.m1 * p.l1^2 / 12;                         % [kg*m^2]
p.J2 = p.m2 * p.l2^2 / 12;                         % [kg*m^2]

% Optional payload. Disabled in the three-rigid-body baseline. When enabled,
% it is counted exactly once at the link-2 tip. Set model to 'point' (J=0)
% or 'rigid' (specified planar COM inertia J).
p.payload.enabled = false;
p.payload.model = 'point';
p.payload.mass = 0.0;        % [kg]
p.payload.J = 0.0;           % [kg*m^2]

validate_parameters(p);
end
