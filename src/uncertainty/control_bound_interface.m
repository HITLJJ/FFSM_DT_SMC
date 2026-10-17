function [bound,cache]=control_bound_interface(mode,z,set,geometry,cache,time,cfg)
%CONTROL_BOUND_INTERFACE Unified nominal torque and uncertainty-bound query.
% Inputs use the Stage-2 reference state z. Sampling-based Taylor values are
% EXPLORATORY and never relabeled as certified. Cache is twin/controller-local.

if nargin<5 || isempty(cache), cache=struct('next_update',-inf,'beta',[0;0]); end
theta=set.nominal(:);
bound.u0=reference_torque_map(z,theta,geometry);
switch lower(mode)
    case 'fixed'
        bound.beta=cfg.fixed_rho(:);
        bound.status='FIXED_ASSUMPTION';
        bound.method='User-configured fixed robust gain.';
    case 'candidate'
        if time+10*eps(max(1,time))>=cache.next_update
            J=parameter_jacobian(z,theta,geometry);
            H=parameter_hessian(z,theta,geometry);
            L=[norm(H(:,:,1),2);norm(H(:,:,2),2)];
            r=norm(set.halfwidth,2);
            cache.beta=abs(J)*set.halfwidth(:)+0.5*L*r^2;
            cache.next_update=time+cfg.bound_update_period;
            cache.last_jacobian=J; cache.last_hessian_norm=L;
        end
        bound.beta=max(cache.beta,0);
        bound.status='EXPLORATORY';
        bound.method='State-dependent local Taylor Hessian candidate (not globally certified).';
    case 'stage2_certified'
        certificate=certified_reference_bound(z,set,geometry);
        bound.beta=certificate.channel_bound;
        bound.status=certificate.status;
        bound.method='Stage-2 global analytic certificate.';
    otherwise
        error('FFSM:UnknownBoundMode','Unknown bound mode %s.',mode);
end
bound.parameter_set=set;
end
