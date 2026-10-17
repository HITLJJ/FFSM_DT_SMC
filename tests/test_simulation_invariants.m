function tests = test_simulation_invariants
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
cfg = stage1_config();
cfg.t_final = 2.0;
cfg.output_dt = 0.02;
testCase.TestData.cfg = cfg;
end

function testExperimentAStatic(testCase)
result = simulate_stage1_case('A',testCase.TestData.cfg);
drift = max(vecnorm((result.x-result.x(1,:)).',2,1));
verifyLessThan(testCase,drift,1e-11);
end

function testExperimentBConservation(testCase)
cfg = testCase.TestData.cfg; result = simulate_stage1_case('B',cfg);
d = result.diagnostics; th = cfg.thresholds;
verifyLessThan(testCase,d.max_linear_momentum_drift,th.momentum);
verifyLessThan(testCase,d.max_angular_momentum_drift,th.momentum);
verifyLessThan(testCase,d.max_energy_drift,th.energy);
verifyLessThan(testCase,d.max_acceleration_error,th.acceleration);
verifyLessThan(testCase,d.max_reconstruction_error,th.reconstruction);
end

function testExperimentCWorkEnergyAndMomentum(testCase)
cfg = testCase.TestData.cfg; result = simulate_stage1_case('C',cfg);
d = result.diagnostics; th = cfg.thresholds;
verifyLessThan(testCase,d.max_linear_momentum_drift,th.momentum);
verifyLessThan(testCase,d.max_angular_momentum_drift,th.momentum);
verifyLessThan(testCase,d.max_work_energy_error,th.work_energy);
verifyGreaterThan(testCase,max(abs(d.T-d.T(1))),1e-8);
end

function testToleranceConvergence(testCase)
cfg = testCase.TestData.cfg; result = simulate_stage1_case('C',cfg);
verifyLessThan(testCase,result.convergence.max_state_difference, ...
    cfg.thresholds.convergence);
end
