function tau = stage1_torque(caseName, t, ~)
%STAGE1_TORQUE Deterministic bounded internal-joint torque profiles [N*m].

switch upper(caseName)
    case {'A','B'}
        tau = zeros(2,1);
    case 'C'
        ramp = 1 - exp(-(t/0.4)^2);
        tau = ramp * [0.35*sin(0.8*t); 0.25*sin(1.1*t + 0.3)];
    otherwise
        error('FFSM:UnknownCase', 'Unknown experiment case %s.', caseName);
end
end
