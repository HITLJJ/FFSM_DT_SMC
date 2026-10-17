function cfg=closed_loop_config()
%CLOSED_LOOP_CONFIG Reproducible Plant-DT-SMC study settings in SI units.
cfg.nominal_parameters=default_parameters();
cfg.initial_set=parameter_set(cfg.nominal_parameters,0.20);
cfg.t_final=6.0; cfg.control_dt=0.01;
cfg.seed=20261017;
cfg.solver='sampled_rk4';
cfg.actuator_limit=[12;12]; % N*m, experimental assumption (not hardware data)
cfg.lambda=[2.5;2.5];
cfg.K=[3.0;3.0];
cfg.phi=[0.04;0.04];
cfg.eta=[0.15;0.15];
cfg.fixed_rho=[1.2;1.2];
cfg.disturbance_bound=[0.30;0.25];
cfg.adaptive.initial_rho=[0.4;0.4];
cfg.adaptive.gamma=[2.0;2.0];
cfg.adaptive.decay=[0.35;0.35];
cfg.adaptive.maximum=[6;6];
cfg.bound_update_period=0.10;
cfg.bound_hessian_mode='local_nominal';
cfg.twin.filter_alpha=1.0;
cfg.twin.update_period=0.10;
cfg.twin.update_method='cumulative_nominal';
cfg.twin.parameter_gain=0.20;
cfg.twin.regularization=0.10;
cfg.twin.information_forgetting=0.995;
cfg.twin.shrink_factor=0.997;
cfg.twin.minimum_relative_halfwidth=0.03;
cfg.twin.maximum_step_fraction=0.01;
root=fileparts(fileparts(mfilename('fullpath')));
cfg.output_dir=fullfile(root,'results','closed_loop');
end
