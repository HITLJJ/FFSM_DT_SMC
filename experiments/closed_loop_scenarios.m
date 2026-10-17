function scenarios=closed_loop_scenarios(cfg)
%CLOSED_LOOP_SCENARIOS Shared Plant conditions for fair controller studies.
theta0=nominal_inertial_parameters(cfg.nominal_parameters);
mismatch=[1.12;0.92;1.18;1.10;0.90;1.15].*theta0;
scenarios(1)=make('nominal',theta0,false);
scenarios(2)=make('parameter_mismatch',mismatch,false);
scenarios(3)=make('mismatch_disturbance',mismatch,true);
end

function s=make(name,theta,disturbed)
s=struct('name',name,'true_theta',theta(:),'disturbed',disturbed);
end
