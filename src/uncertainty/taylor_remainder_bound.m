function out = taylor_remainder_bound(z, set, geometry, sampleCount, seed)
%TAYLOR_REMAINDER_BOUND Support-function--Taylor CANDIDATE bound.
% Hessian norms are sampled numeric references, not global upper envelopes;
% therefore this function can return CANDIDATE or FAILED, never CERTIFIED.

if nargin<4, sampleCount=12; end
if nargin<5, seed=20261017; end
tStart=tic;
J=parameter_jacobian(z,set.nominal,geometry);
linear=support_function_bound(J,set);
rng(seed,'twister');
points=[set.nominal, set.lower, set.upper];
points=[points, set.lower+(set.upper-set.lower).*rand(6,sampleCount)];
L=zeros(2,1);
for k=1:size(points,2)
    H=parameter_hessian(z,points(:,k),geometry);
    for channel=1:2
        L(channel)=max(L(channel),norm(H(:,:,channel),2));
    end
end
r=norm(set.halfwidth,2);
out.parameter_jacobian=J;
out.linear_radius=linear.linear_radius;
out.sampled_hessian_norm=L;
out.parameter_radius=r;
out.channel_bound=out.linear_radius+0.5*L*r^2;
out.status='CANDIDATE';
out.reason='Hessian maximum is sampled and has no validated global enclosure.';
out.sample_count=size(points,2);
out.seed=seed;
out.elapsed_seconds=toc(tStart);
end
