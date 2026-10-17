function create_stage1_animation(result, p, outputFile)
%CREATE_STAGE1_ANIMATION Save planar robot motion as Motion JPEG AVI.

t = result.t; x = result.x;
n = numel(t); stride = max(1,round(n/300));
allPoints = zeros(2,3,n);
for k = 1:n
    kin = planar_kinematics(x(k,1:5).',p);
    allPoints(:,:,k) = [kin.shoulder kin.elbow kin.end_effector];
end
xlimv = [min(allPoints(1,:),[],'all') max(allPoints(1,:),[],'all')] + [-0.75 0.75];
ylimv = [min(allPoints(2,:),[],'all') max(allPoints(2,:),[],'all')] + [-0.75 0.75];

writer = VideoWriter(outputFile,'Motion JPEG AVI');
writer.FrameRate = 30; writer.Quality = 95; open(writer);
fig = figure('Visible','off','Color','w');
cleanup = onCleanup(@() close_resources(writer,fig));
for k = 1:stride:n
    Q = x(k,1:5).'; kin = planar_kinematics(Q,p);
    cla; hold on; grid on; axis equal;
    psi = Q(3); R = [cos(psi) -sin(psi);sin(psi) cos(psi)];
    corners = R*[p.base_half_length p.base_half_length -p.base_half_length -p.base_half_length; ...
        p.base_half_width -p.base_half_width -p.base_half_width p.base_half_width] + Q(1:2);
    patch(corners(1,:),corners(2,:),[0.90 0.70 0.15],'EdgeColor',[0.2 0.2 0.2]);
    plot([kin.shoulder(1) kin.elbow(1) kin.end_effector(1)], ...
        [kin.shoulder(2) kin.elbow(2) kin.end_effector(2)],'-o','LineWidth',2);
    plot(squeeze(allPoints(1,3,1:k)),squeeze(allPoints(2,3,1:k)),':k');
    xlim(xlimv); ylim(ylimv); xlabel('X [m]'); ylabel('Y [m]');
    title(sprintf('Experiment %s, t = %.2f s',result.case_name,t(k)));
    drawnow; writeVideo(writer,getframe(fig));
end
close(writer); close(fig); clear cleanup;
end

function close_resources(writer,fig)
try, close(writer); catch, end
try, close(fig); catch, end
end
