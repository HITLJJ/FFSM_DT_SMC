function summary=run_stage3a_hessian_certification()
%RUN_STAGE3A_HESSIAN_CERTIFICATION Compare certified and exploratory bounds.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src'))); addpath(fullfile(root,'config'));
cfg=stage3a_config(); if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

row=0; records=struct([]); failures={};
for stateIndex=1:numel(cfg.states)
    z=cfg.states(stateIndex);
    coefficientTimer=tic; coeff=fixed_state_block_norms(z.q,cfg.geometry);
    coefficientSeconds=toc(coefficientTimer);
    for radiusIndex=1:numel(cfg.relative_radii)
        radius=cfg.relative_radii(radiusIndex); set=parameter_set(cfg.geometry,radius);
        baseline=certified_reference_bound(z,set,cfg.geometry);
        candidate=taylor_remainder_bound(z,set,cfg.geometry,cfg.hessian_samples, ...
            cfg.seed+100*stateIndex+radiusIndex);
        strict=certified_hessian_bound(z,set,cfg.geometry,coeff);
        empirical=empirical_bound_search(z,set,cfg.geometry,cfg.random_samples, ...
            cfg.seed+1000*stateIndex+radiusIndex);
        if ~strcmp(strict.status,'CERTIFIED')
            failures{end+1}=struct('state',z.name,'radius',radius, ...
                'reason',strict.failure_reason); %#ok<AGROW>
        end
        for channel=1:2
            row=row+1;
            records(row).state=z.name;
            records(row).relative_radius=radius;
            records(row).channel=channel;
            records(row).status=strict.status;
            records(row).stage2_certified=baseline.channel_bound(channel);
            records(row).stage2_candidate=candidate.channel_bound(channel);
            records(row).stage3_certified=value_or_nan(strict,'bound',channel);
            records(row).stage3_elementwise=value_or_nan(strict,'elementwise_remainder',channel)+ ...
                value_or_nan(strict,'linear_part',channel);
            records(row).stage3_raw_spectral=value_or_nan(strict,'raw_spectral_remainder',channel)+ ...
                value_or_nan(strict,'linear_part',channel);
            records(row).stage3_normalized_spectral= ...
                value_or_nan(strict,'normalized_spectral_remainder',channel)+ ...
                value_or_nan(strict,'linear_part',channel);
            records(row).linear_part=value_or_nan(strict,'linear_part',channel);
            records(row).remainder_part=value_or_nan(strict,'remainder_part',channel);
            if isfield(strict,'chosen_remainder_method')
                records(row).chosen_remainder_method=strict.chosen_remainder_method{channel};
            else
                records(row).chosen_remainder_method='';
            end
            records(row).remainder_fraction=records(row).remainder_part/ ...
                records(row).stage3_certified;
            records(row).empirical_max=empirical.max_abs_error(channel);
            records(row).stage2_ratio=baseline.channel_bound(channel)/empirical.max_abs_error(channel);
            records(row).stage3_ratio=records(row).stage3_certified/empirical.max_abs_error(channel);
            records(row).reduction_factor=baseline.channel_bound(channel)/ ...
                records(row).stage3_certified;
            records(row).coefficient_seconds=coefficientSeconds;
            records(row).certificate_seconds=value_or_nan(strict,'runtime',1);
            records(row).candidate_seconds=candidate.elapsed_seconds;
            records(row).empirical_seconds=empirical.elapsed_seconds;
            records(row).failure_reason=getfield_or(strict,'failure_reason','');
        end
    end
end
t=struct2table(records); writetable(t,fullfile(cfg.output_dir,'stage3a_summary.csv'));

stress=run_stress(cfg); writetable(stress,fullfile(cfg.output_dir,'stress_summary.csv'));
save(fullfile(cfg.output_dir,'stage3a_raw.mat'),'cfg','records','stress','failures');
make_plots(t,stress,cfg.output_dir);
summary=struct('channel_cases',height(t),'certified_cases',sum(strcmp(t.status,'CERTIFIED')), ...
    'failures',numel(failures),'median_reduction',median(t.reduction_factor,'omitnan'), ...
    'max_reduction',max(t.reduction_factor,[],'omitnan'), ...
    'median_stage3_ratio',median(t.stage3_ratio,'omitnan'), ...
    'max_stage3_ratio',max(t.stage3_ratio,[],'omitnan'), ...
    'stress_certified',sum(strcmp(stress.status,'CERTIFIED')), ...
    'stress_total',height(stress),'output_dir',cfg.output_dir);
write_summary(summary,t,cfg.output_dir);
end

function stress=run_stress(cfg)
rows=struct([]); set=parameter_set(cfg.geometry,0.20);
for n=1:numel(cfg.stress)
    z=cfg.stress(n); coeff=fixed_state_block_norms(z.q,cfg.geometry);
    c=certified_hessian_bound(z,set,cfg.geometry,coeff);
    e=empirical_bound_search(z,set,cfg.geometry,80,cfg.seed+9000+n);
    rows(n).state=z.name;
    rows(n).status=c.status;
    rows(n).bound_channel_1=value_or_nan(c,'bound',1);
    rows(n).bound_channel_2=value_or_nan(c,'bound',2);
    rows(n).empirical_channel_1=e.max_abs_error(1);
    rows(n).empirical_channel_2=e.max_abs_error(2);
    rows(n).runtime=value_or_nan(c,'runtime',1);
    rows(n).failure_reason=getfield_or(c,'failure_reason','');
end
stress=struct2table(rows);
end

function make_plots(t,stress,out)
x=(1:height(t)).'; labels=strcat(string(t.state),"-",string(100*t.relative_radius),"%-c",string(t.channel));
f=figure('Visible','off','Color','w','Position',[80,80,1250,520]);
semilogy(x,t.empirical_max,'ko'); hold on; semilogy(x,t.stage2_certified,'r-');
semilogy(x,t.stage2_candidate,'b--'); semilogy(x,t.stage3_certified,'g-','LineWidth',1.3); grid on;
legend('EMPIRICAL','Stage 2 CERTIFIED','Stage 2 CANDIDATE','Stage 3A CERTIFIED','Location','best');
title('Full reference-torque uncertainty bounds'); ylabel('N m'); xticks(x); xticklabels(labels); xtickangle(75);
export_plot(f,out,'all_bound_comparison'); close(f);

f=figure('Visible','off','Color','w');
boxchart(categorical(t.relative_radius),t.stage3_certified); set(gca,'YScale','log'); grid on;
xlabel('Relative box half width'); ylabel('Stage 3A CERTIFIED bound [N m]');
title('Certified width versus uncertainty'); export_plot(f,out,'width_vs_uncertainty'); close(f);

f=figure('Visible','off','Color','w');
semilogy(x,t.stage2_certified./t.empirical_max,'r-'); hold on;
semilogy(x,t.stage3_certified./t.empirical_max,'g-');
semilogy(x,t.stage2_candidate./t.empirical_max,'b--'); grid on;
legend('Stage 2 CERTIFIED','Stage 3A CERTIFIED','Stage 2 CANDIDATE');
ylabel('Bound / empirical maximum'); xlabel('Case'); title('Conservatism comparison');
export_plot(f,out,'conservatism_comparison'); close(f);

f=figure('Visible','off','Color','w');
methodValues=[median(t.stage3_elementwise),median(t.stage3_raw_spectral), ...
    median(t.stage3_normalized_spectral)];
bar(methodValues); set(gca,'YScale','log');
ylim([min(methodValues)/2,max(methodValues)*2]); grid on;
set(gca,'XTickLabel',{'Elementwise','Raw spectral','Normalized spectral'});
ylabel('Median complete bound [N m]'); title('Strict Hessian remainder algorithms');
export_plot(f,out,'remainder_methods'); close(f);

f=figure('Visible','off','Color','w');
bar([median(t.coefficient_seconds),median(t.certificate_seconds), ...
    median(t.empirical_seconds)]); grid on;
set(gca,'XTickLabel',{'Trig enclosure','Hessian certificate','Empirical search'});
ylabel('Median time [s]'); title('Computation time'); export_plot(f,out,'computation_time'); close(f);

f=figure('Visible','off','Color','w');
bar([stress.bound_channel_1,stress.bound_channel_2]); set(gca,'YScale','log'); grid on;
xticks(1:height(stress)); xticklabels(stress.state); xtickangle(20);
set(gca,'TickLabelInterpreter','none');
ylabel('CERTIFIED bound [N m]'); title('Stage 3A stress cases');
legend('channel 1','channel 2'); export_plot(f,out,'stress_cases'); close(f);
end

function write_summary(s,t,out)
fid=fopen(fullfile(out,'summary.txt'),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'Stage 3A Hessian certification\n');
fprintf(fid,'Certified channel cases: %d/%d\n',s.certified_cases,s.channel_cases);
fprintf(fid,'Certification failures: %d\n',s.failures);
fprintf(fid,'Stage 2 / Stage 3A reduction median/max: %.9g / %.9g\n',s.median_reduction,s.max_reduction);
fprintf(fid,'Stage 3A / empirical ratio median/max: %.9g / %.9g\n',s.median_stage3_ratio,s.max_stage3_ratio);
fprintf(fid,'Remainder fraction median/max: %.9g / %.9g\n',median(t.remainder_fraction),max(t.remainder_fraction));
fprintf(fid,'Stress certified: %d/%d\n',s.stress_certified,s.stress_total);
fprintf(fid,'Maximum certificate-only time: %.9g s\n',max(t.certificate_seconds));
fprintf(fid,'Maximum fixed-state coefficient time: %.9g s\n',max(t.coefficient_seconds));
fprintf(fid,'Sampling and fmincon remain empirical counterexample searches only.\n');
end

function value=value_or_nan(s,field,index)
if isfield(s,field) && ~isempty(s.(field)), v=s.(field); value=v(index); else, value=NaN; end
end
function value=getfield_or(s,field,fallback)
if isfield(s,field), value=s.(field); else, value=fallback; end
end
function export_plot(f,out,name)
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',180);
savefig(f,fullfile(out,[name '.fig']));
end
