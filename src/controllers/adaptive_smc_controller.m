function [tau,state,details]=adaptive_smc_controller(signal,bound,cfg,state,dt)
%ADAPTIVE_SMC_CONTROLLER Explicit adaptive-gain SMC baseline wrapper.
state.gain_mode='adaptive';
[tau,state,details]=smc_controller(signal,bound,cfg,state,dt);
details.adaptation_status='EXPLORATORY_BASELINE';
end
