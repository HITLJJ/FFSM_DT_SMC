function twin=twin_initialize(geometry,initialSet,options,mode)
%TWIN_INITIALIZE Create a state-isolated digital twin.
% The twin receives copied measurements only. It cannot write Plant state or
% access a hidden true parameter. Its outputs may influence the controller.

twin.geometry=geometry;
twin.theta_hat=initialSet.nominal(:);
twin.parameter_set=initialSet;
twin.initial_set=initialSet;
twin.options=options;
twin.mode=mode;
twin.filtered_q=[]; twin.filtered_qdot=[];
twin.previous_q=[]; twin.previous_qdot=[]; twin.previous_time=[];
twin.latest_raw_measurement=[];
twin.last_update_time=-inf;
twin.update_count=0;
twin.version=1;
twin.validity_status='INITIALIZED';
twin.last_residual=[NaN;NaN];
twin.last_parameter_step=zeros(6,1);
twin.information_matrix=zeros(6);
twin.information_vector=zeros(6,1);
end
