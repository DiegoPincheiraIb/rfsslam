close all

%%hfig = figure('Renderer','OpenGL');
%%set(gcf, 'Color', 'w');
%%title('Groundturth robot trajectory and landmark positions');
%%plot(gt_pose(1,:), gt_pose(2,:), 'r--');
%%hold on
%%plot(zeros(1, length(gt_lmk(1,:))), gt_lmk(1,:), 'k.');
%axis square
%%grid on
%%set(gca, 'XLim', get(gca, 'XLim') + [-1, 1] );

hfig = figure('Renderer','OpenGL');
set(gcf, 'Color', 'w');
plot(gt_lmk(1,:), zeros(1, length(gt_lmk(1,:))), 'k.');
hold on
h_robotPos = plot(gt_pose(2,1), 0, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
h_drPos = plot(dr_pose(2,1), -0.1, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
h_particlePos = plot(x_i(1, :, 1), 0, 'm.');
grid on
set(gca, 'XLim', get(gca, 'XLim') + [-1, 1] );
set(gca, 'YLim', [-1, 15] );
xLim = get(gca, 'XLim');
vEvalPts = xLim(1) : 0.005 : xLim(2);
v = vEvalPts * 0;
h_map = plot(vEvalPts, v, 'b-');

pause(0.1);
meas_idx = 1;

for k = 1 : length(gt_pose)
    
    disp(k);

    delete( h_robotPos )
    delete( h_drPos )
    delete( h_particlePos );
    delete(findobj('Color','b'))
    
    x = gt_pose(2,k);
    h_robotPos = plot(x, 0, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    
    dr_x = dr_pose(2,k);
    h_drPos = plot(dr_x, -0.1, 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
    
    h_particlePos = plot(x_i(1, :, k), 0, 'm.');

    line([x+3.1 x+3.1], [-1.5 0], 'Color', 'b');
    line([x+2.9 x+2.9], [-1.5 0], 'Color', 'b');
    line([x+0.4 x+0.4], [-1.5 0], 'Color', 'b');
    line([x+0.6 x+0.6], [-1.5 0], 'Color', 'b');
    line([x-3.1 x-3.1], [-1.5 0], 'Color', 'b');
    line([x-2.9 x-2.9], [-1.5 0], 'Color', 'b');
    line([x-0.4 x-0.4], [-1.5 0], 'Color', 'b');
    line([x-0.6 x-0.6], [-1.5 0], 'Color', 'b');
    if(meas_idx <= length(meas))
        while( meas(1, meas_idx) == k )
            r = meas(2, meas_idx);
            line([x+r x+r], [-0.5 0.5], 'Color', 'b');
            meas_idx = meas_idx + 1;
            if(meas_idx > length(meas))
                break
            end
        end
    end
    
    highest_weight_i = 0;
    highest_weight = 0;
    for i = 1:nParticles{1}
        if x_i(2, i, k) > highest_weight
            highest_weight = x_i(2, i, k);
            highest_weight_i = i;
        end
    end
    
    v = vEvalPts * 0;
    for i = 1:nParticles{1}
        [nLm tmp] = size(landmarkEst{i,k});
        for m = 1:nLm
            u = landmarkEst{i,k}{m,1};
            cov = landmarkEst{i,k}{m,2};
            w = landmarkEst{i,k}{m,3};
            v = v + pdf('norm', vEvalPts, u, cov*3) * w;
            %if( w > 0.5)
            %    line([u, u],[0 15], 'Color', 'b')
            %end
        end         
    end
    v = v / double(nParticles{1});
    h_map = plot(vEvalPts, v, 'b-');
    
    %export_fig(sprintf('results/anim/%06d.png',k), hfig);
    pause(0.05)
    
end