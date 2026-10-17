function tests=test_smc_closed_loop
tests=functiontests(localfunctions);
end

function setupOnce(~)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src'))); addpath(fullfile(root,'config'));
addpath(fullfile(root,'experiments'));
end

function testNominalFullPlantTracking(testCase)
cfg=closed_loop_config(); cfg.t_final=1.0;
s=closed_loop_scenarios(cfg); c=closed_loop_controllers();
r=simulate_closed_loop_case(s(1),c(1),cfg);
verifyEqual(testCase,r.status,'COMPLETED');
verifyLessThan(testCase,r.metrics.position_rmse,1e-4);
verifyLessThanOrEqual(testCase,max(abs(r.tau),[],'all'),max(cfg.actuator_limit));
verifyTrue(testCase,all(isfinite(r.x),'all'));
end

function testNoBaseActuatorAndMomentum(testCase)
cfg=closed_loop_config(); cfg.t_final=0.5;
s=closed_loop_scenarios(cfg); c=closed_loop_controllers();
r=simulate_closed_loop_case(s(2),c(1),cfg);
p=inertial_parameters_to_model(s(2).true_theta,cfg.nominal_parameters);
worst=0;
for k=1:size(r.x,1)
    me=momentum_energy(r.x(k,1:5).',r.x(k,6:10).',p);
    worst=max(worst,norm([me.P;me.LO]));
end
verifyLessThan(testCase,worst,1e-8);
end
