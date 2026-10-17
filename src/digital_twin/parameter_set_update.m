function [twin,info]=parameter_set_update(twin,z,measuredTorque)
%PARAMETER_SET_UPDATE Projected normalized Gauss-Newton exploratory update.
% Uses only synchronized state/rate changes and applied joint torque. Because
% acceleration differentiation, disturbance and noise bounds are not certified,
% set containment is EXPLORATORY and no true-parameter inclusion is claimed.

theta=twin.theta_hat; theta0=twin.initial_set.nominal(:); geometry=twin.geometry;
scale=theta0;
switch lower(twin.options.update_method)
    case 'single_sample'
        predicted=reference_torque_map(z,theta,geometry);
        residual=measuredTorque(:)-predicted;
        J=parameter_jacobian(z,theta,geometry); Js=J*diag(scale);
        deltaXi=twin.options.parameter_gain*((Js.'*Js+ ...
            twin.options.regularization*eye(6))\(Js.'*residual));
    case 'cumulative_nominal'
        predicted=reference_torque_map(z,theta0,geometry);
        residual=measuredTorque(:)-predicted;
        J=parameter_jacobian(z,theta0,geometry); Js=J*diag(scale);
        forget=twin.options.information_forgetting;
        twin.information_matrix=forget*twin.information_matrix+Js.'*Js;
        twin.information_vector=forget*twin.information_vector+Js.'*residual;
        targetXi=(twin.information_matrix+twin.options.regularization*eye(6))\ ...
            twin.information_vector;
        radius=(twin.initial_set.upper-twin.initial_set.nominal)./scale;
        targetXi=max(-radius,min(radius,targetXi));
        desired=theta0+scale.*targetXi;
        deltaXi=twin.options.parameter_gain*(desired-theta)./scale;
    otherwise
        error('FFSM:UnknownTwinUpdate','Unknown twin update method %s.', ...
            twin.options.update_method);
end
maxStep=twin.options.maximum_step_fraction*ones(6,1);
% parameter_gain is applied once in each update law above. Applying it again
% here would silently square the configured gain.
deltaXi=max(-maxStep,min(maxStep,deltaXi));
candidate=theta+scale.*deltaXi;
candidate=max(twin.initial_set.lower,min(twin.initial_set.upper,candidate));

oldHalf=twin.parameter_set.halfwidth(:);
minimum=twin.options.minimum_relative_halfwidth*scale;
newHalf=max(minimum,twin.options.shrink_factor*oldHalf);
lowerBound=max(twin.initial_set.lower,candidate-newHalf);
upperBound=min(twin.initial_set.upper,candidate+newHalf);
candidate=min(upperBound,max(lowerBound,candidate));
twin.theta_hat=candidate;
twin.parameter_set.nominal=candidate;
twin.parameter_set.lower=lowerBound;
twin.parameter_set.upper=upperBound;
twin.parameter_set.halfwidth=max(candidate-lowerBound,upperBound-candidate);
twin.parameter_set.description='EXPLORATORY online projected parameter box';
twin.parameter_set.is_hardware_confidence_set=false;
twin.parameter_set.certification_status='EXPLORATORY';
info=struct('residual',residual,'parameter_step',candidate-theta, ...
    'jacobian',J,'status','EXPLORATORY');
end
