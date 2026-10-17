function [tau,state,details]=smc_controller(signal,bound,cfg,state,dt)
%SMC_CONTROLLER First-order saturated SMC for two FFSM joints.
% signal fields q, qdot, qd, qd_dot [rad,rad/s]. bound.u0 and beta are N*m.
% The output is actuator-limited torque [N*m]. No base actuation is created.

e=signal.q(:)-signal.qd(:);
edot=signal.qdot(:)-signal.qd_dot(:);
s=edot+cfg.lambda(:).*e;
switch lower(state.gain_mode)
    case 'fixed'
        rho=cfg.fixed_rho(:);
    case 'adaptive'
        drho=cfg.adaptive.gamma(:).*abs(s)-cfg.adaptive.decay(:).* ...
            (state.rho-cfg.adaptive.initial_rho(:));
        state.rho=min(cfg.adaptive.maximum(:),max(cfg.adaptive.initial_rho(:), ...
            state.rho+dt*drho));
        rho=state.rho;
    case 'bound'
        rho=bound.beta(:)+cfg.disturbance_bound(:)+cfg.eta(:);
    otherwise
        error('FFSM:UnknownGainMode','Unknown SMC gain mode %s.',state.gain_mode);
end
switching=max(-1,min(1,s./cfg.phi(:)));
tauRaw=bound.u0(:)-cfg.K(:).*s-rho.*switching;
limit=cfg.actuator_limit(:);
tau=max(-limit,min(limit,tauRaw));
details=struct('e',e,'edot',edot,'s',s,'rho',rho,'switching',switching, ...
    'tau_raw',tauRaw,'saturated',abs(tauRaw)>limit,'bound_status',bound.status);
end
