function m=compute_closed_loop_metrics(result,cfg)
%COMPUTE_CLOSED_LOOP_METRICS Reproducible tracking, effort and twin metrics.
valid=all(isfinite(result.error),2); assert(any(valid),'FFSM:NoValidSamples','No valid samples.');
e=result.error(valid,:);
tau=result.tau(valid,:); s=result.s(valid,:); beta=result.beta(valid,:);
modelError=abs(result.model_error(valid,:)); cycle=result.cycle_time(valid);
m.position_rmse_channels=sqrt(mean(e.^2,1));
m.position_rmse=sqrt(mean(e.^2,'all'));
qdot=result.x(1:size(result.qdot_desired,1),9:10); qdot=qdot(valid,:);
m.velocity_rmse=sqrt(mean((qdot-result.qdot_desired(valid,:)).^2,'all'));
base=result.x(1:size(result.error,1),1:3); base=base(valid,:);
m.base_translation_peak=max(vecnorm((base(:,1:2)-base(1,1:2)).',2,1));
m.base_yaw_peak=max(abs(base(:,3)-base(1,3)));
m.torque_rms=sqrt(mean(tau.^2,'all'));
m.torque_peak=max(abs(tau),[],'all');
m.torque_total_variation=sum(abs(diff(tau,1,1)),'all');
m.saturation_ratio=mean(result.saturated(valid,:),'all');
m.sliding_rms=sqrt(mean(s.^2,'all'));
m.final_error_norm=norm(e(end,:));
m.model_error_rmse=sqrt(mean(modelError.^2,'all'));
positive=modelError>1e-9;
ratios=beta(positive)./modelError(positive);
if isempty(ratios)
    m.bound_error_ratio_median=NaN;
    m.bound_error_ratio_max=NaN;
    m.bound_violation_count=0;
else
    m.bound_error_ratio_median=median(ratios);
    m.bound_error_ratio_max=max(1./ratios);
    m.bound_violation_count=sum(ratios<1);
end
m.twin_update_count=result.twin_final.update_count;
m.parameter_error_initial=norm(cfg.initial_set.nominal-result.scenario.true_theta);
m.parameter_error_final=norm(result.twin_final.theta_hat-result.scenario.true_theta);
m.parameter_set_width_initial=sum(2*cfg.initial_set.halfwidth);
m.parameter_set_width_final=sum(result.twin_final.parameter_set.upper- ...
    result.twin_final.parameter_set.lower);
m.cycle_time_mean=mean(cycle); m.cycle_time_max=max(cycle);
m.failed=~strcmp(result.status,'COMPLETED');
m.final_velocity_snapshot=result.x(find(valid,1,'last'),9:10);
end
