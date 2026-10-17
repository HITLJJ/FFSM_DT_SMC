function cfg = stage1_config()
%STAGE1_CONFIG Numerical settings for stage-1 dynamics validation.
% Units: SI (m, kg, s, rad, N*m).

cfg.parameters = default_parameters();
cfg.t_final = 10.0;
cfg.output_dt = 0.01;
cfg.solver = 'ode113';

cfg.strict.RelTol = 1e-10;
cfg.strict.AbsTol = 1e-12;
cfg.strict.MaxStep = 0.01;

cfg.coarse.RelTol = 1e-7;
cfg.coarse.AbsTol = 1e-9;
cfg.coarse.MaxStep = 0.05;

cfg.thresholds.symmetry = 1e-12;
cfg.thresholds.positive_definite = 1e-8;
cfg.thresholds.acceleration = 1e-8;
cfg.thresholds.momentum = 1e-7;
cfg.thresholds.energy = 1e-7;
cfg.thresholds.work_energy = 1e-6;
cfg.thresholds.reconstruction = 1e-8;
cfg.thresholds.finite_difference = 1e-6;
cfg.thresholds.convergence = 1e-5;

cfg.output_dir = fullfile(project_root(), 'results', 'stage1');
end

function root = project_root()
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
end
