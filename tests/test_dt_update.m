function tests=test_dt_update
tests=functiontests(localfunctions);
end

function setupOnce(~)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src'))); addpath(fullfile(root,'config'));
addpath(fullfile(root,'experiments'));
end

function testRawMeasurementNotMutated(testCase)
cfg=closed_loop_config(); packet=plant_measurement_packet(0,(1:10).',1);
original=packet; twin=twin_initialize(cfg.nominal_parameters,cfg.initial_set,cfg.twin,'fixed_model');
[updatedTwin,output]=twin_step(twin,packet,[NaN;NaN]);
verifyEqual(testCase,packet,original);
verifyEqual(testCase,output.q_estimate,original.raw_Q(4:5));
verifyFalse(testCase,isfield(updatedTwin,'true_theta'));
end

function testTwinChangesControlThroughInterface(testCase)
cfg=closed_loop_config(); cfg.t_final=3.0;
s=closed_loop_scenarios(cfg); c=closed_loop_controllers();
fixed=simulate_closed_loop_case(s(2),c(3),cfg);
updated=simulate_closed_loop_case(s(2),c(4),cfg);
verifyGreaterThan(testCase,updated.twin_final.update_count,0);
verifyGreaterThan(testCase,norm(updated.u0-fixed.u0,'fro'),1e-6);
verifyGreaterThan(testCase,norm(updated.tau-fixed.tau,'fro'),1e-7);
verifyEqual(testCase,updated.twin_final.parameter_set.certification_status,'EXPLORATORY');
end
