function runs=run_closed_loop_dt_smc(cfg)
%RUN_CLOSED_LOOP_DT_SMC Experiment 4: isolate online-twin benefit without disturbance.
if nargin<1, cfg=closed_loop_config(); end
s=closed_loop_scenarios(cfg); c=closed_loop_controllers();
runs=cell(2,1);
runs{1}=simulate_closed_loop_case(s(2),c(3),cfg); runs{1}.run_name='E4_C3_fixed_twin';
runs{2}=simulate_closed_loop_case(s(2),c(4),cfg); runs{2}.run_name='E4_C4_online_twin';
end
