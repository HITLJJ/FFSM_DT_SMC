function out = empirical_bound_search(z, set, geometry, sampleCount, seed)
%EMPIRICAL_BOUND_SEARCH Search for large fixed-z torque-map deviations.
% Corners, pseudorandom samples, and (when available) local fmincon searches
% are evaluated. Results are empirical lower bounds on the worst case only;
% they are never mathematical certificates. theta units are kg and kg*m^2,
% while the returned channel deviations are N*m.

if nargin<4, sampleCount=256; end
if nargin<5, seed=20261017; end
validateattributes(sampleCount,{'numeric'},{'scalar','integer','nonnegative'});
tStart=tic;
rng(seed,'twister');

bits=dec2bin(0:63,6).'-'0';
corners=set.lower+(set.upper-set.lower).*bits;
randomPoints=set.lower+(set.upper-set.lower).*rand(6,sampleCount);
points=[set.nominal,corners,randomPoints];
F0=reference_torque_map(z,set.nominal,geometry);
[maxAbs,argmax,source]=scan_points(points,F0,z,geometry,'grid');

usedFmincon=false;
fminconRuns=0;
if exist('fmincon','file')==2 && license('test','Optimization_Toolbox')
    usedFmincon=true;
    options=optimoptions('fmincon','Display','none','Algorithm','sqp', ...
        'MaxFunctionEvaluations',500,'OptimalityTolerance',1e-10);
    starts=[set.nominal,argmax,randomPoints(:,1:min(2,sampleCount))];
    for channel=1:2
        for signValue=[-1,1]
            for k=1:size(starts,2)
                [thetaStar,~,exitflag]=fmincon(@(theta) channel_value( ...
                    theta,z,geometry,F0,channel,signValue),starts(:,k),[],[],[],[], ...
                    set.lower,set.upper,[],options);
                fminconRuns=fminconRuns+1;
                if exitflag>=0
                    [newMax,newArg,newSource]=scan_points(thetaStar,F0,z,geometry,'fmincon');
                    improve=newMax>maxAbs;
                    maxAbs(improve)=newMax(improve);
                    argmax(:,improve)=newArg(:,improve);
                    source(improve)=newSource(improve);
                end
            end
        end
    end
end

out.status='EMPIRICAL';
out.max_abs_error=maxAbs;
out.argmax=argmax;
out.argmax_source=source;
out.nominal_value=F0;
out.sample_count=size(points,2);
out.random_count=sampleCount;
out.corner_count=64;
out.seed=seed;
out.used_fmincon=usedFmincon;
out.fmincon_runs=fminconRuns;
out.elapsed_seconds=toc(tStart);
out.warning=['Observed maxima are lower bounds on the true worst case; ', ...
    'sampling and local optimization do not prove validity.'];
end

function value=channel_value(theta,z,geometry,F0,channel,signValue)
difference=reference_torque_map(z,theta,geometry)-F0;
value=-signValue*difference(channel);
end

function [maxAbs,argmax,source]=scan_points(points,F0,z,geometry,label)
errors=zeros(2,size(points,2));
for k=1:size(points,2)
    errors(:,k)=abs(reference_torque_map(z,points(:,k),geometry)-F0);
end
[maxAbs,index]=max(errors,[],2);
argmax=[points(:,index(1)),points(:,index(2))];
source=repmat({label},2,1);
end
