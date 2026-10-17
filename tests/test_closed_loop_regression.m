function tests=test_closed_loop_regression
tests=functiontests(localfunctions);
end

function setupOnce(~)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src'))); addpath(fullfile(root,'config'));
addpath(fullfile(root,'experiments'));
end

function testAllControllersCompleteCommonScenario(testCase)
cfg=closed_loop_config(); cfg.t_final=1.0;
s=closed_loop_scenarios(cfg); c=closed_loop_controllers();
for k=1:numel(c)
    r=simulate_closed_loop_case(s(3),c(k),cfg);
    verifyEqual(testCase,r.status,'COMPLETED');
    verifyFalse(testCase,r.metrics.failed);
    verifyLessThanOrEqual(testCase,r.metrics.torque_peak,max(cfg.actuator_limit));
    verifyLessThan(testCase,r.metrics.position_rmse,0.05);
end
end

function testCandidateStatusRemainsExploratory(testCase)
cfg=closed_loop_config(); cfg.t_final=0.2;
s=closed_loop_scenarios(cfg); c=closed_loop_controllers();
r=simulate_closed_loop_case(s(2),c(3),cfg);
verifyTrue(testCase,all(strcmp(r.bound_status,'EXPLORATORY')));
end
