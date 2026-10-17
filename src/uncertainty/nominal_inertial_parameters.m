function theta = nominal_inertial_parameters(p)
%NOMINAL_INERTIAL_PARAMETERS Return [m0;m1;m2;J0;J1;J2] in SI units.
% Masses use kg and planar COM inertias use kg*m^2.

theta = [p.m0;p.m1;p.m2;p.J0;p.J1;p.J2];
end
