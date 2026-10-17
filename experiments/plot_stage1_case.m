function plot_stage1_case(result, outputDir)
%PLOT_STAGE1_CASE Save native MATLAB figures and PNG validation plots.

name = ['experiment_' lower(result.case_name)];
t = result.t;
x = result.x;
d = result.diagnostics;

fig = figure('Visible','off','Color','w');
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
nexttile(tl); plot(t,x(:,1:2),'LineWidth',1.2); grid on;
xlabel('Time [s]'); ylabel('Base position [m]'); legend('x_b','y_b','Location','best');
nexttile(tl); plot(t,x(:,3),'LineWidth',1.2); grid on;
xlabel('Time [s]'); ylabel('\psi [rad]');
nexttile(tl); plot(t,x(:,4:5),'LineWidth',1.2); grid on;
xlabel('Time [s]'); ylabel('Joint angle [rad]'); legend('q_1','q_2','Location','best');
nexttile(tl); plot(t,x(:,9:10),'LineWidth',1.2); grid on;
xlabel('Time [s]'); ylabel('Joint rate [rad/s]'); legend('dq_1','dq_2','Location','best');
title(tl,['Experiment ' result.case_name ': states']);
save_both(fig,outputDir,[name '_states']); close(fig);

fig = figure('Visible','off','Color','w');
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
nexttile(tl); plot(t,d.delta_P,'LineWidth',1.2); grid on;
xlabel('Time [s]'); ylabel('\DeltaP [kg m/s]'); legend('\DeltaP_x','\DeltaP_y','Location','best');
nexttile(tl); plot(t,d.delta_LO,'LineWidth',1.2); grid on;
xlabel('Time [s]'); ylabel('\DeltaL_O [kg m^2/s]');
title(tl,['Experiment ' result.case_name ': momentum error']);
save_both(fig,outputDir,[name '_momentum']); close(fig);

fig = figure('Visible','off','Color','w');
plot(t,d.T-d.T(1),'LineWidth',1.2); hold on;
plot(t,d.work,'--','LineWidth',1.2); plot(t,d.work_energy_error,':','LineWidth',1.2); grid on;
xlabel('Time [s]'); ylabel('Energy [J]');
legend('\DeltaT','\int dq^T\tau dt','work-energy error','Location','best');
title(['Experiment ' result.case_name ': work-energy consistency']);
save_both(fig,outputDir,[name '_energy']); close(fig);

fig = figure('Visible','off','Color','w');
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
nexttile(tl); semilogy(t,max(d.acceleration_error,eps),'LineWidth',1.2); grid on;
ylabel('||ddq_f-ddq_r||_2 [rad/s^2]');
nexttile(tl); semilogy(t,max(d.base_velocity_reconstruction_error,eps),'LineWidth',1.2); grid on;
ylabel('||db-db_r||_2');
nexttile(tl); plot(t,d.min_eig_H,'LineWidth',1.2); hold on;
plot(t,d.min_eig_Mr,'--','LineWidth',1.2); grid on;
xlabel('Time [s]'); ylabel('Minimum eigenvalue'); legend('H','M_r','Location','best');
title(tl,['Experiment ' result.case_name ': model consistency']);
save_both(fig,outputDir,[name '_consistency']); close(fig);
end

function save_both(fig,outputDir,baseName)
if ~exist(outputDir,'dir'), mkdir(outputDir); end
savefig(fig,fullfile(outputDir,[baseName '.fig']));
exportgraphics(fig,fullfile(outputDir,[baseName '.png']),'Resolution',200);
end
