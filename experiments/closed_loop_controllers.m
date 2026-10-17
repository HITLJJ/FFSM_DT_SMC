function controllers=closed_loop_controllers()
%CLOSED_LOOP_CONTROLLERS Controller definitions used in every fair comparison.
controllers(1)=make('C1_fixed_smc','fixed','fixed','fixed_model');
controllers(2)=make('C2_adaptive_smc','fixed','adaptive','fixed_model');
controllers(3)=make('C3_candidate_smc','candidate','bound','fixed_model');
controllers(4)=make('C4_dt_candidate_smc','candidate','bound','online_update');
end

function c=make(name,boundMode,gainMode,twinMode)
c=struct('name',name,'bound_mode',boundMode,'gain_mode',gainMode,'twin_mode',twinMode);
end
