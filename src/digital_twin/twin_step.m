function [twin,output]=twin_step(twin,measurement,previousAppliedTorque)
%TWIN_STEP Causal measurement synchronization and optional parameter update.
% measurement.raw_Q/raw_Qdot are immutable Plant observations. Filtering,
% estimation and prediction occur on twin-owned copies. The twin never writes
% back to measurement or Plant; influence is through output -> controller.

assert(isstruct(measurement) && isfield(measurement,'raw_Q'), ...
    'FFSM:InvalidMeasurement','A timestamped raw Plant measurement is required.');
alpha=twin.options.filter_alpha;
q=measurement.raw_Q(4:5); qdot=measurement.raw_Qdot(4:5);
if isempty(twin.filtered_q)
    twin.filtered_q=q; twin.filtered_qdot=qdot;
else
    twin.filtered_q=alpha*q+(1-alpha)*twin.filtered_q;
    twin.filtered_qdot=alpha*qdot+(1-alpha)*twin.filtered_qdot;
end
twin.latest_raw_measurement=measurement;

canUpdate=strcmpi(twin.mode,'online_update') && ~isempty(twin.previous_time) && ...
    measurement.time-twin.last_update_time>=twin.options.update_period-1e-12 && ...
    all(isfinite(previousAppliedTorque));
if canUpdate
    dt=measurement.time-twin.previous_time;
    if dt>0
        qdd=(twin.filtered_qdot-twin.previous_qdot)/dt;
        qmid=0.5*(twin.filtered_q+twin.previous_q);
        qdotmid=0.5*(twin.filtered_qdot+twin.previous_qdot);
        z=struct('q',qmid,'qdot',qdotmid, ...
            'qr_dot',qdotmid,'qr_ddot',qdd);
        [twin,update]=parameter_set_update(twin,z,previousAppliedTorque(:));
        twin.last_residual=update.residual;
        twin.last_parameter_step=update.parameter_step;
        twin.last_update_time=measurement.time;
        twin.update_count=twin.update_count+1;
        twin.version=twin.version+1;
    end
end
twin.previous_q=twin.filtered_q;
twin.previous_qdot=twin.filtered_qdot;
twin.previous_time=measurement.time;
twin.validity_status=twin_validity_check(twin,measurement.time);
output=struct('q_estimate',twin.filtered_q,'qdot_estimate',twin.filtered_qdot, ...
    'theta_hat',twin.theta_hat,'parameter_set',twin.parameter_set, ...
    'validity_status',twin.validity_status,'version',twin.version, ...
    'update_count',twin.update_count,'last_residual',twin.last_residual);
end
