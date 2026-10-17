function result=simulate_closed_loop_case(scenario,controller,cfg)
%SIMULATE_CLOSED_LOOP_CASE Sampled-data full 5-DOF Plant-DT-SMC simulation.
% The Plant uses hidden scenario.true_theta. Twin/controller receive only raw
% measurement packets and applied commands. All coordinates and outputs are SI.

rng(cfg.seed,'twister');
pTrue=inertial_parameters_to_model(scenario.true_theta,cfg.nominal_parameters);
ref0=closed_loop_reference(0); Q0=[0;0;0;ref0.q];
Qdot0=[zero_momentum_base_velocity(Q0,ref0.qdot,pTrue);ref0.qdot];
x=[Q0;Qdot0];
t=(0:cfg.control_dt:cfg.t_final).'; n=numel(t); steps=n-1;
twin=twin_initialize(cfg.nominal_parameters,cfg.initial_set,cfg.twin,controller.twin_mode);
controllerState=struct('gain_mode',controller.gain_mode,'rho',cfg.adaptive.initial_rho);
boundCache=[]; previousTau=[NaN;NaN];

result.t=t; result.x=nan(n,10); result.x(1,:)=x.';
fields={'q_desired','qdot_desired','tau','tau_raw','disturbance','error','s','rho', ...
    'beta','u0','model_error','theta_hat','set_halfwidth','cycle_time'};
widths=[2,2,2,2,2,2,2,2,2,2,2,6,6,1];
for k=1:numel(fields), result.(fields{k})=nan(steps,widths(k)); end
result.saturated=false(steps,2); result.bound_status=cell(steps,1);
result.twin_validity=cell(steps,1); result.twin_version=zeros(steps,1);
result.raw_measurement_q=nan(steps,2); result.twin_q_estimate=nan(steps,2);
status='COMPLETED'; failureReason='';
for k=1:steps
    cycleTimer=tic;
    measurement=plant_measurement_packet(t(k),x,k);
    [twin,twinOut]=twin_step(twin,measurement,previousTau);
    ref=closed_loop_reference(t(k));
    q=twinOut.q_estimate; qdot=twinOut.qdot_estimate;
    e=q-ref.q; edot=qdot-ref.qdot;
    z=struct('q',q,'qdot',qdot, ...
        'qr_dot',ref.qdot-cfg.lambda(:).*e, ...
        'qr_ddot',ref.qddot-cfg.lambda(:).*edot);
    [bound,boundCache]=control_bound_interface(controller.bound_mode,z, ...
        twinOut.parameter_set,cfg.nominal_parameters,boundCache,t(k),cfg);
    signal=struct('q',q,'qdot',qdot,'qd',ref.q,'qd_dot',ref.qdot);
    if strcmp(controller.gain_mode,'adaptive')
        [tau,controllerState,control]=adaptive_smc_controller( ...
            signal,bound,cfg,controllerState,cfg.control_dt);
    else
        [tau,controllerState,control]=smc_controller( ...
            signal,bound,cfg,controllerState,cfg.control_dt);
    end
    disturbance=matched_disturbance(t(k),scenario);
    actualF=reference_torque_map(z,scenario.true_theta,cfg.nominal_parameters);

    result.q_desired(k,:)=ref.q.'; result.qdot_desired(k,:)=ref.qdot.';
    result.tau(k,:)=tau.'; result.tau_raw(k,:)=control.tau_raw.';
    result.disturbance(k,:)=disturbance.'; result.error(k,:)=control.e.';
    result.s(k,:)=control.s.'; result.rho(k,:)=control.rho.';
    result.beta(k,:)=bound.beta.'; result.u0(k,:)=bound.u0.';
    result.model_error(k,:)=(actualF-bound.u0).';
    result.theta_hat(k,:)=twinOut.theta_hat.';
    result.set_halfwidth(k,:)=twinOut.parameter_set.halfwidth.';
    result.saturated(k,:)=control.saturated.';
    result.bound_status{k}=bound.status; result.twin_validity{k}=twinOut.validity_status;
    result.twin_version(k)=twinOut.version;
    result.raw_measurement_q(k,:)=measurement.raw_Q(4:5).';
    result.twin_q_estimate(k,:)=twinOut.q_estimate.';

    x=rk4_step(t(k),x,tau,scenario,pTrue,cfg.control_dt);
    result.cycle_time(k)=toc(cycleTimer);
    if any(~isfinite(x)) || norm(x,inf)>1e4
        status='FAILED'; failureReason=sprintf('Nonfinite or divergent state at step %d.',k);
        result.x(k+1,:)=x.'; break
    end
    result.x(k+1,:)=x.';
    previousTau=tau;
end
result.status=status; result.failure_reason=failureReason;
result.scenario=scenario; result.controller=controller;
result.twin_final=twin; result.true_parameters_for_evaluation=scenario.true_theta;
result.metrics=compute_closed_loop_metrics(result,cfg);
end

function next=rk4_step(time,x,tau,scenario,p,dt)
f=@(tt,xx) plant_rhs(tt,xx,tau,scenario,p);
k1=f(time,x); k2=f(time+dt/2,x+dt*k1/2);
k3=f(time+dt/2,x+dt*k2/2); k4=f(time+dt,x+dt*k3);
next=x+dt*(k1+2*k2+2*k3+k4)/6;
end

function dx=plant_rhs(time,x,tau,scenario,p)
d=matched_disturbance(time,scenario);
Qdd=full_acceleration(x(1:5),x(6:10),tau(:)+d,p);
dx=[x(6:10);Qdd];
end

function d=matched_disturbance(time,scenario)
if scenario.disturbed
    d=[0.25*sin(1.7*time)+0.05*sin(4.3*time);0.20*cos(1.3*time)];
else
    d=[0;0];
end
end
