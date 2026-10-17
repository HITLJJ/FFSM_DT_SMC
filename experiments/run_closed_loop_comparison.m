function summary=run_closed_loop_comparison()
%RUN_CLOSED_LOOP_COMPARISON End-to-end Plant-DT-SMC experiments and iterations.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src'))); addpath(fullfile(root,'config'));
cfg=closed_loop_config(); if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
s=closed_loop_scenarios(cfg); c=closed_loop_controllers();

runs=[run_closed_loop_baselines(cfg);run_closed_loop_dt_smc(cfg)];
hardC4=simulate_closed_loop_case(s(3),c(4),cfg); hardC4.run_name='E3_C4_dt_disturbance';
runs{end+1}=hardC4;
metrics=metrics_table(runs);
writetable(metrics,fullfile(cfg.output_dir,'closed_loop_metrics.csv'));

[iterationRuns,iterations]=run_iterations(cfg,s,c);
writetable(iterations,fullfile(cfg.output_dir,'iteration_history.csv'));
save(fullfile(cfg.output_dir,'closed_loop_raw.mat'),'cfg','s','c','runs', ...
    'metrics','iterationRuns','iterations','-v7.3');
make_closed_loop_plots(runs,metrics,iterations,cfg.output_dir);

summary=struct(); summary.run_count=numel(runs);
summary.completed=sum(strcmp(metrics.status,'COMPLETED'));
summary.failed=summary.run_count-summary.completed;
summary.best_rmse=min(metrics.position_rmse);
summary.best_controller=metrics.controller{find(metrics.position_rmse==summary.best_rmse,1)};
summary.dt_changes_u0=norm(runs{8}.u0-runs{7}.u0,'fro');
summary.dt_changes_beta=norm(runs{8}.beta-runs{7}.beta,'fro');
summary.dt_changes_tau=norm(runs{8}.tau-runs{7}.tau,'fro');
summary.output_dir=cfg.output_dir;
write_text_summary(summary,metrics,iterations,cfg.output_dir);
end

function [allRuns,tableOut]=run_iterations(cfg,s,c)
definitions={ ...
    'round1_single_sample','single_sample',0.03,0.001, ...
        'Single-sample Gauss-Newton; exposed disturbance attribution failure.'; ...
    'round2_cumulative','cumulative_nominal',0.20,0.10, ...
        'Cumulative nominal linearization; retained for disturbance-free update benefit.'; ...
    'round3_regularized','cumulative_nominal',0.10,1.00, ...
        'Stronger regularization; reduced drift but largely removed independent benefit.'};
allRuns=cell(size(definitions,1),2); row=0; records=struct([]);
for round=1:size(definitions,1)
    local=cfg; local.twin.update_method=definitions{round,2};
    local.twin.parameter_gain=definitions{round,3};
    local.twin.regularization=definitions{round,4};
    for scenarioIndex=2:3
        r=simulate_closed_loop_case(s(scenarioIndex),c(4),local);
        allRuns{round,scenarioIndex-1}=r; row=row+1;
        records(row).round=definitions{round,1};
        records(row).scenario=s(scenarioIndex).name;
        records(row).method=definitions{round,2};
        records(row).gain=definitions{round,3};
        records(row).regularization=definitions{round,4};
        records(row).reason=definitions{round,5};
        records(row).position_rmse=r.metrics.position_rmse;
        records(row).model_error_rmse=r.metrics.model_error_rmse;
        records(row).torque_rms=r.metrics.torque_rms;
        records(row).torque_total_variation=r.metrics.torque_total_variation;
        records(row).parameter_error_initial=r.metrics.parameter_error_initial;
        records(row).parameter_error_final=r.metrics.parameter_error_final;
        records(row).update_count=r.metrics.twin_update_count;
        records(row).status=r.status;
    end
end
tableOut=struct2table(records);
end

function t=metrics_table(runs)
records=struct([]);
for k=1:numel(runs)
    r=runs{k}; m=r.metrics;
    records(k).run=r.run_name; records(k).scenario=r.scenario.name;
    records(k).controller=r.controller.name; records(k).status=r.status;
    records(k).position_rmse=m.position_rmse;
    records(k).velocity_rmse=m.velocity_rmse;
    records(k).base_translation_peak=m.base_translation_peak;
    records(k).base_yaw_peak=m.base_yaw_peak;
    records(k).torque_rms=m.torque_rms; records(k).torque_peak=m.torque_peak;
    records(k).torque_total_variation=m.torque_total_variation;
    records(k).saturation_ratio=m.saturation_ratio;
    records(k).sliding_rms=m.sliding_rms;
    records(k).model_error_rmse=m.model_error_rmse;
    records(k).bound_error_ratio_median=m.bound_error_ratio_median;
    records(k).bound_error_ratio_max=m.bound_error_ratio_max;
    records(k).bound_violation_count=m.bound_violation_count;
    records(k).twin_update_count=m.twin_update_count;
    records(k).parameter_error_initial=m.parameter_error_initial;
    records(k).parameter_error_final=m.parameter_error_final;
    records(k).parameter_set_width_initial=m.parameter_set_width_initial;
    records(k).parameter_set_width_final=m.parameter_set_width_final;
    records(k).cycle_time_mean=m.cycle_time_mean;
    records(k).cycle_time_max=m.cycle_time_max;
end
t=struct2table(records);
end

function make_closed_loop_plots(runs,metrics,iterations,out)
hard={runs{4},runs{5},runs{6},runs{9}};
f=figure('Visible','off','Color','w','Position',[80,80,1100,700]);
tiledlayout(2,1);
for joint=1:2
    nexttile; plot(hard{1}.t(1:end-1),hard{1}.q_desired(:,joint),'k--','LineWidth',1.2); hold on;
    for k=1:4, plot(hard{k}.t,hard{k}.x(:,3+joint),'LineWidth',1.0); end
    grid on; ylabel(sprintf('q_%d [rad]',joint));
    if joint==1, title('Matched-disturbance tracking on the full free-floating Plant'); end
end
xlabel('Time [s]'); legend('reference','C1','C2','C3','C4','Location','best');
export_plot(f,out,'tracking_comparison'); close(f);

f=figure('Visible','off','Color','w','Position',[80,80,1100,700]); tiledlayout(2,2);
for k=1:4
    nexttile; plot(hard{k}.t(1:end-1),hard{k}.tau); grid on;
    yline(12,'r:'); yline(-12,'r:'); title(hard{k}.controller.name,'Interpreter','none');
    xlabel('Time [s]'); ylabel('Torque [N m]');
end
export_plot(f,out,'torque_comparison'); close(f);

mask=strcmp(metrics.scenario,'mismatch_disturbance'); common=metrics(mask,:);
f=figure('Visible','off','Color','w','Position',[80,80,1000,520]);
values=[common.position_rmse,common.torque_rms,common.torque_total_variation/20, ...
    common.sliding_rms]; bar(values); grid on;
xticklabels(common.controller); set(gca,'TickLabelInterpreter','none');
legend('position RMSE','torque RMS','TV / 20','sliding RMS','Location','best');
title('Fair common-scenario metrics (different physical units shown separately)');
export_plot(f,out,'common_metrics'); close(f);

fixed=runs{7}; updated=runs{8};
f=figure('Visible','off','Color','w','Position',[80,80,1100,700]); tiledlayout(2,1);
nexttile; plot(updated.t(1:end-1),vecnorm(updated.model_error,2,2)); hold on;
plot(fixed.t(1:end-1),vecnorm(fixed.model_error,2,2),'--'); grid on;
ylabel('||F(true)-u_0|| [N m]'); legend('online twin','fixed twin'); title('Digital-twin causal influence');
nexttile; plot(updated.t(1:end-1),updated.theta_hat); grid on;
ylabel('Parameter estimate'); xlabel('Time [s]');
export_plot(f,out,'digital_twin_effect'); close(f);

f=figure('Visible','off','Color','w','Position',[80,80,1000,550]);
actual=abs(updated.model_error); plot(updated.t(1:end-1),actual,'LineWidth',1); hold on;
plot(updated.t(1:end-1),updated.beta,'--','LineWidth',1); grid on;
xlabel('Time [s]'); ylabel('Torque error / candidate bound [N m]');
legend('|Delta F_1|','|Delta F_2|','candidate beta_1','candidate beta_2');
title('EXPLORATORY candidate bound audit'); export_plot(f,out,'candidate_bound_audit'); close(f);

f=figure('Visible','off','Color','w','Position',[80,80,1000,520]);
groups=categorical(strcat(iterations.round," / ",iterations.scenario));
bar(groups,[iterations.parameter_error_final,iterations.model_error_rmse]); grid on;
set(gca,'TickLabelInterpreter','none'); xtickangle(25);
legend('final parameter error','model torque RMSE'); title('Three recorded twin-update iterations');
export_plot(f,out,'iteration_history'); close(f);
end

function write_text_summary(s,m,it,out)
fid=fopen(fullfile(out,'summary.txt'),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'Closed-loop Plant-DT-SMC research summary\n');
fprintf(fid,'Completed runs: %d/%d; failures: %d\n',s.completed,s.run_count,s.failed);
fprintf(fid,'Best recorded position RMSE: %.9g (%s)\n',s.best_rmse,s.best_controller);
fprintf(fid,'DT causal changes ||du0||/||dbeta||/||dtau||: %.9g / %.9g / %.9g\n', ...
    s.dt_changes_u0,s.dt_changes_beta,s.dt_changes_tau);
fprintf(fid,'Common disturbed comparison:\n');
mask=strcmp(m.scenario,'mismatch_disturbance');
for k=find(mask).', fprintf(fid,'  %s RMSE %.9g torqueRMS %.9g TV %.9g\n', ...
        m.controller{k},m.position_rmse(k),m.torque_rms(k),m.torque_total_variation(k)); end
fprintf(fid,'Iteration outcomes (including adverse results):\n');
for k=1:height(it)
    fprintf(fid,'  %s %s pErr %.9g->%.9g modelRMSE %.9g\n', ...
        it.round{k},it.scenario{k},it.parameter_error_initial(k), ...
        it.parameter_error_final(k),it.model_error_rmse(k));
end
fprintf(fid,'Candidate-bound and online-update conclusions are EXPLORATORY.\n');
end

function export_plot(f,out,name)
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',180);
savefig(f,fullfile(out,[name '.fig']));
end
