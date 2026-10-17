function summary = run_stage1_validation()
%RUN_STAGE1_VALIDATION Run experiments A-C and produce reproducible artifacts.
% Entry point: from the repository root, run
%   addpath(genpath(pwd)); summary = run_stage1_validation;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
addpath(fullfile(root,'config'));
addpath(fullfile(root,'experiments'));
cfg = stage1_config();
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

cases = {'A','B','C'};
summary = struct();
for i = 1:numel(cases)
    result = simulate_stage1_case(cases{i},cfg);
    save(fullfile(cfg.output_dir,['experiment_' lower(cases{i}) '.mat']), ...
        'result','cfg','-v7.3');
    plot_stage1_case(result,cfg.output_dir);
    summary.(cases{i}) = scalar_summary(result);
    if strcmp(cases{i},'C')
        create_stage1_animation(result,cfg.parameters, ...
            fullfile(cfg.output_dir,'experiment_c_animation.avi'));
    end
end
save(fullfile(cfg.output_dir,'stage1_summary.mat'),'summary','cfg');
write_summary_text(summary,cfg,fullfile(cfg.output_dir,'stage1_summary.txt'));
disp(summary);
end

function s = scalar_summary(r)
d = r.diagnostics;
s.max_linear_momentum_drift = d.max_linear_momentum_drift;
s.max_angular_momentum_drift = d.max_angular_momentum_drift;
s.max_energy_drift = d.max_energy_drift;
s.max_work_energy_error = d.max_work_energy_error;
s.max_acceleration_error = d.max_acceleration_error;
s.max_reconstruction_error = d.max_reconstruction_error;
s.max_symmetry_error = d.max_symmetry_error;
s.minimum_eigenvalue_H = d.minimum_eigenvalue_H;
s.minimum_eigenvalue_Mr = d.minimum_eigenvalue_Mr;
s.max_state_convergence_difference = r.convergence.max_state_difference;
end

function write_summary_text(summary,cfg,fileName)
fid = fopen(fileName,'w');
assert(fid >= 0,'FFSM:IO','Could not create summary file.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'Stage 1 dynamics validation\n');
fprintf(fid,'Solver: %s\nStrict: RelTol=%g AbsTol=%g MaxStep=%g\n', ...
    cfg.solver,cfg.strict.RelTol,cfg.strict.AbsTol,cfg.strict.MaxStep);
names = fieldnames(summary);
for i = 1:numel(names)
    fprintf(fid,'\nExperiment %s\n',names{i});
    values = summary.(names{i}); fields = fieldnames(values);
    for j = 1:numel(fields)
        fprintf(fid,'%s = %.16g\n',fields{j},values.(fields{j}));
    end
end
end
