function set = parameter_set(geometry, relativeRadius)
%PARAMETER_SET Axis-aligned Stage-2 inertial-parameter test box.
% relativeRadius is a scalar fraction (0.05, 0.10 or 0.20). The returned box
% is an algorithm-test domain, not a hardware confidence interval and not a
% complete rigid-body physical-consistency set.

validateattributes(relativeRadius,{'numeric'},{'real','finite','scalar','>',0,'<',1});
theta0 = nominal_inertial_parameters(geometry);
set.kind = 'axis_aligned_box';
set.relative_radius = relativeRadius;
set.nominal = theta0;
set.halfwidth = relativeRadius*theta0;
set.lower = theta0-set.halfwidth;
set.upper = theta0+set.halfwidth;
set.parameter_names = {'m0','m1','m2','J0','J1','J2'};
set.units = {'kg','kg','kg','kg*m^2','kg*m^2','kg*m^2'};
set.is_hardware_confidence_set = false;
set.includes_joint_physical_constraints = false;
set.description = sprintf('Independent nominal +/- %.0f%% inertial test box',100*relativeRadius);
assert(all(set.lower>0),'FFSM:InvalidParameterSet','Parameter box violates positivity.');
end
