function cfg=stage3a_config()
%STAGE3A_CONFIG Reproducible Hessian-certification and stress-test settings.
cfg=stage2_config();
cfg.seed=20261018;
cfg.output_dir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','stage3a');
cfg.stress(1)=state('box_boundary',[0.45;-1.1],[0.6;-0.5],[0.5;-0.4],[1.0;-0.8]);
cfg.stress(2)=state('near_singular_configuration',[pi-1e-8;-pi+2e-8], ...
    [0.7;-0.6],[0.5;-0.45],[1.2;-1.0]);
cfg.stress(3)=state('large_reference_velocity',[-1.3;1.8],[4;-3],[5;-4],[1;-0.8]);
cfg.stress(4)=state('large_reference_acceleration',[0.9;-1.7],[0.5;0.4], ...
    [0.4;0.3],[10;-8]);
end

function z=state(name,q,qdot,qrDot,qrDDot)
z=struct('name',name,'q',q(:),'qdot',qdot(:), ...
    'qr_dot',qrDot(:),'qr_ddot',qrDDot(:));
end
