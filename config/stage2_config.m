function cfg = stage2_config()
%STAGE2_CONFIG Fixed states and numerical settings for uncertainty validation.
% All states use joint coordinates q [rad], physical joint rates qdot [rad/s],
% reference rates qr_dot [rad/s], and reference accelerations qr_ddot [rad/s^2].

cfg.geometry=default_parameters();
cfg.relative_radii=[0.05,0.10,0.20];
cfg.random_samples=180;
cfg.hessian_samples=5;
cfg.seed=20261017;

cfg.states(1)=make_state('extended',[0;0],[0.12;-0.08],[0.10;-0.06],[0.20;-0.15]);
cfg.states(2)=make_state('folded',[pi;0],[-0.15;0.12],[-0.10;0.08],[-0.25;0.18]);
cfg.states(3)=make_state('bent',[0.7;-1.25],[0.30;-0.22],[0.24;-0.17],[0.45;-0.35]);
cfg.states(4)=make_state('high_velocity',[-0.9;1.4],[1.8;-1.5],[1.5;-1.2],[0.35;-0.30]);
cfg.states(5)=make_state('high_acceleration',[1.2;-0.6],[0.25;0.18],[0.20;0.15],[3.0;-2.5]);
cfg.states(6)=make_state('strong_base_reaction',[1.55;-2.2],[1.3;1.1],[-1.0;1.4],[2.2;2.6]);

root=fileparts(fileparts(mfilename('fullpath')));
cfg.output_dir=fullfile(root,'results','stage2');
end

function z=make_state(name,q,qdot,qrDot,qrDDot)
z=struct('name',name,'q',q(:),'qdot',qdot(:), ...
    'qr_dot',qrDot(:),'qr_ddot',qrDDot(:));
end
