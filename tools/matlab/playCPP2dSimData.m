close all

hfig = figure('Renderer','OpenGL');
set(gcf, 'Color', 'w');
title('Groundturth robot trajectory and landmark positions');
plot(gt_pose(2,:), gt_pose(3,:), 'r--');
hold on
plot(gt_lmk(1,:), gt_lmk(2,:), 'k.');
axis equal
grid on
set(gca, 'XLim', get(gca, 'XLim') + [-1, 1] );
set(gca, 'YLim', get(gca, 'YLim') + [-1, 1] );
   
h_robotPos = plot(gt_pose(2,1), gt_pose(3,1), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
h_robotHdg = line([gt_pose(2,1) gt_pose(2,1)+0.5*cos(gt_pose(4,1))], [gt_pose(3,1) gt_pose(3,1)+0.5*sin(gt_pose(4,1))], 'Color', 'k');
h_drPos = plot(dr_pose(2,1), dr_pose(3,1), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
h_drHdg = line([dr_pose(2,1) dr_pose(2,1)+0.5*cos(dr_pose(4,1))], [dr_pose(3,1) dr_pose(3,1)+0.5*sin(dr_pose(4,1))], 'Color', 'k');
h_particlePos = plot(x_i(1, :, 1), x_i(2, :, 1), 'm.');

pause(0.005);
meas_idx = 1;

for k = 1 : length(gt_pose)

    delete( h_robotPos )
    delete( h_robotHdg )
    delete( h_drPos )
    delete( h_drHdg )
    delete( h_particlePos );
    delete(findobj('Color','b'))
    
    t = gt_pose(1,k);
    x = gt_pose(2,k);
    y = gt_pose(3,k);
    z = gt_pose(4,k);
    h_robotPos = plot(x, y, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    h_robotHdg = line([x x+0.5*cos(z)], [y y+0.5*sin(z)], 'Color', 'k');
    
    dr_x = dr_pose(2,k);
    dr_y = dr_pose(3,k);
    dr_z = dr_pose(4,k);
    h_drPos = plot(dr_x, dr_y, 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
    h_drHdg = line([dr_x dr_x+0.5*cos(dr_z)], [dr_y dr_y+0.5*sin(dr_z)], 'Color', 'k');
    
    h_particlePos = plot(x_i(1, :, k), x_i(2, :, k), 'm.');
    
    if(meas_idx <= length(meas))
        while( meas(1, meas_idx) == t )
            r = meas(2, meas_idx);
            b = meas(3, meas_idx);
            mx = x + r*cos(b + z);
            my = y + r*sin(b + z);
            line([x mx], [y my], 'Color', 'b');
            meas_idx = meas_idx + 1;
            if(meas_idx > length(meas))
                break
            end
        end
    end
    
    highest_weight_i = 0;
    highest_weight = 0;
    for i = 1:nParticles{1}
        if x_i(4, i, k) > highest_weight
            highest_weight = x_i(4, i, k);
            highest_weight_i = i;
        end
    end
    for i = highest_weight_i
        [nLm tmp] = size(landmarkEst{i,k});
        for m = 1:nLm
            u = landmarkEst{i,k}{m,1};
            cov = landmarkEst{i,k}{m,2};
            w = landmarkEst{i,k}{m,3};
            [evec, eval] = eig(cov);
            axes_length = 3 * sqrt(diag(eval));
            if(axes_length(2) > axes_length(1))
                angle = atan2(evec(2,1), evec(1,1));
            else
                angle = atan2(evec(2,2), evec(1,2));
            end
            h_ellipse = ellipse(axes_length(1), axes_length(2), angle, u(1), u(2));
        end
    end

    
    %export_fig(sprintf('results/anim/%06d.png',k), hfig);
    pause(0.1)
    
end