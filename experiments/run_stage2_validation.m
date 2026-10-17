function summary=run_stage2_validation()
%RUN_STAGE2_VALIDATION Validate strict and candidate Stage-2 uncertainty bounds.
% Generated data are empirical diagnostics unless explicitly marked CERTIFIED.
% No controller, estimator, or digital-twin update is performed here.

root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
addpath(fullfile(root,'config'));
cfg=stage2_config();
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

row=0; records=struct([]); counterexamples={}; sensitivities=zeros(2,6,numel(cfg.states));
for stateIndex=1:numel(cfg.states)
    z=cfg.states(stateIndex);
    sensitivities(:,:,stateIndex)=parameter_jacobian( ...
        z,nominal_inertial_parameters(cfg.geometry),cfg.geometry);
    for radiusIndex=1:numel(cfg.relative_radii)
        radius=cfg.relative_radii(radiusIndex);
        set=parameter_set(cfg.geometry,radius);
        strictTimer=tic;
        certified=certified_reference_bound(z,set,cfg.geometry);
        certifiedTime=toc(strictTimer);
        candidate=taylor_remainder_bound(z,set,cfg.geometry, ...
            cfg.hessian_samples,cfg.seed+100*stateIndex+radiusIndex);
        empirical=empirical_bound_search(z,set,cfg.geometry, ...
            cfg.random_samples,cfg.seed+1000*stateIndex+radiusIndex);

        candidateStatus=repmat({candidate.status},2,1);
        failed=empirical.max_abs_error>candidate.channel_bound.*(1+1e-10);
        candidateStatus(failed)={'FAILED'};
        certifiedFailed=empirical.max_abs_error>certified.channel_bound.*(1+1e-12);
        assert(~any(certifiedFailed),'FFSM:CertifiedBoundViolation', ...
            'An empirical point exceeded a CERTIFIED bound.');
        for channel=1:2
            row=row+1;
            records(row).state=z.name;
            records(row).relative_radius=radius;
            records(row).channel=channel;
            records(row).empirical_max=empirical.max_abs_error(channel);
            records(row).certified_bound=certified.channel_bound(channel);
            records(row).candidate_bound=candidate.channel_bound(channel);
            records(row).candidate_status=candidateStatus{channel};
            records(row).certified_status=certified.status;
            records(row).inertia_bound=certified.inertia_component;
            records(row).coriolis_bound=certified.coriolis_component;
            records(row).certified_ratio=certified.channel_bound(channel)/ ...
                max(empirical.max_abs_error(channel),realmin);
            records(row).candidate_ratio=candidate.channel_bound(channel)/ ...
                max(empirical.max_abs_error(channel),realmin);
            records(row).certified_seconds=certifiedTime;
            records(row).candidate_seconds=candidate.elapsed_seconds;
            records(row).empirical_seconds=empirical.elapsed_seconds;
            records(row).fmincon_used=empirical.used_fmincon;
            if failed(channel)
                counterexamples{end+1}=struct('state',z.name,'radius',radius, ...
                    'channel',channel,'theta',empirical.argmax(:,channel), ...
                    'observed_error',empirical.max_abs_error(channel), ...
                    'candidate_bound',candidate.channel_bound(channel)); %#ok<AGROW>
            end
        end
    end
end

tableOut=struct2table(records);
writetable(tableOut,fullfile(cfg.output_dir,'bound_summary.csv'));
save(fullfile(cfg.output_dir,'stage2_raw.mat'),'cfg','records','counterexamples', ...
    'sensitivities');
make_plots(tableOut,sensitivities,cfg);

summary=struct();
summary.total_cases=height(tableOut);
summary.certified_failures=sum(tableOut.empirical_max>tableOut.certified_bound);
summary.candidate_failures=sum(strcmp(tableOut.candidate_status,'FAILED'));
summary.max_certified_ratio=max(tableOut.certified_ratio);
summary.median_certified_ratio=median(tableOut.certified_ratio);
summary.max_candidate_ratio=max(tableOut.candidate_ratio);
summary.median_candidate_ratio=median(tableOut.candidate_ratio);
summary.fmincon_used=all(tableOut.fmincon_used);
summary.output_dir=cfg.output_dir;
fid=fopen(fullfile(cfg.output_dir,'summary.txt'),'w');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'Stage 2 validation\n');
fprintf(fid,'Channel cases: %d\n',summary.total_cases);
fprintf(fid,'Certified empirical violations: %d\n',summary.certified_failures);
fprintf(fid,'Taylor candidate failures: %d\n',summary.candidate_failures);
fprintf(fid,'Certified conservatism ratio median/max: %.9g / %.9g\n', ...
    summary.median_certified_ratio,summary.max_certified_ratio);
fprintf(fid,'Candidate ratio median/max: %.9g / %.9g\n', ...
    summary.median_candidate_ratio,summary.max_candidate_ratio);
fprintf(fid,'fmincon used: %d\n',summary.fmincon_used);
fprintf(fid,['Interpretation: sampling/local optimization are empirical only; ', ...
    'CERTIFIED labels derive exclusively from analytic exact-rational envelopes.\n']);
end

function make_plots(t,sensitivities,cfg)
out=cfg.output_dir; x=(1:height(t)).';
labels=strcat(string(t.state),"-",string(100*t.relative_radius),"%-c",string(t.channel));

f=figure('Visible','off','Color','w','Position',[100,100,1200,500]);
semilogy(x,t.empirical_max,'ko','MarkerSize',4); hold on;
semilogy(x,t.certified_bound,'r-','LineWidth',1.2);
semilogy(x,t.candidate_bound,'b--','LineWidth',1.2); grid on;
xlabel('Fixed-state / box / channel case'); ylabel('Torque deviation or bound [N m]');
title('Observed errors and uncertainty bounds (CERTIFIED versus CANDIDATE)');
legend('EMPIRICAL observed max','CERTIFIED analytic bound','Taylor CANDIDATE','Location','best');
xticks(x); xticklabels(labels); xtickangle(75);
export_plot(f,out,'bound_comparison'); close(f);

f=figure('Visible','off','Color','w','Position',[100,100,900,500]);
for k=1:numel(cfg.relative_radii)
    mask=t.relative_radius==cfg.relative_radii(k);
    values=t.certified_bound(mask);
    plot(k*ones(size(values)),values,'ro','MarkerFaceColor',[0.85 0.2 0.2]); hold on;
end
set(gca,'XTick',1:numel(cfg.relative_radii), ...
    'XTickLabel',compose('%g%%',100*cfg.relative_radii)); grid on;
xlabel('Independent parameter-box half width'); ylabel('CERTIFIED bound [N m]');
title('Certified full-map bound growth with uncertainty width');
export_plot(f,out,'bound_vs_uncertainty'); close(f);

f=figure('Visible','off','Color','w','Position',[100,100,950,500]);
ratio=t.certified_bound./max(t.candidate_bound,realmin);
semilogy(x,ratio,'m.-','MarkerSize',10); yline(1,'k:'); grid on;
xlabel('Case'); ylabel('CERTIFIED / Taylor candidate width');
title('Conservatism comparison (Taylor remains CANDIDATE or FAILED)');
export_plot(f,out,'certified_vs_candidate'); close(f);

f=figure('Visible','off','Color','w','Position',[100,100,900,500]);
imagesc(mean(abs(sensitivities),3)); colorbar;
xticks(1:6); xticklabels({'m0','m1','m2','J0','J1','J2'});
yticks(1:2); yticklabels({'F1','F2'});
xlabel('Inertial parameter'); ylabel('Torque-map channel');
title('Mean absolute nominal sensitivity (NUMERICAL)');
export_plot(f,out,'parameter_sensitivity'); close(f);

f=figure('Visible','off','Color','w','Position',[100,100,850,500]);
bar([median(t.certified_seconds),median(t.candidate_seconds),median(t.empirical_seconds)]);
set(gca,'XTickLabel',{'CERTIFIED analytic','Taylor CANDIDATE','EMPIRICAL search'}); grid on;
ylabel('Median wall time [s]'); title('Computation-time comparison');
export_plot(f,out,'computation_time'); close(f);

f=figure('Visible','off','Color','w','Position',[100,100,700,450]);
counts=zeros(numel(cfg.relative_radii),2);
for k=1:numel(cfg.relative_radii)
    mask=t.relative_radius==cfg.relative_radii(k);
    counts(k,1)=sum(t.empirical_max(mask)>t.certified_bound(mask));
    counts(k,2)=sum(strcmp(t.candidate_status(mask),'FAILED'));
end
bar(counts); grid on;
set(gca,'XTickLabel',compose('%g%%',100*cfg.relative_radii));
xlabel('Box half width'); ylabel('Number of empirical violations');
title('Violation audit: certification versus candidate');
legend('CERTIFIED','Taylor FAILED','Location','best');
export_plot(f,out,'violation_audit'); close(f);
end

function export_plot(f,out,name)
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',180);
savefig(f,fullfile(out,[name '.fig']));
end
