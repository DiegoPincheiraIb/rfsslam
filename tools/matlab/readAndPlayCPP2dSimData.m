% Read the data files generated by c++ program sim2d 
% located in [project home]/data/
% Filenames: gtPose.dat, gtLandmark.dat, odometry.dat, measurement.dat
% Keith Leung 2013

addpath '~/src/Matlab/export_fig'

clear all
close all

%logDir = '../../data/rbphdslam/';
%logDir = '../../data/fastslam/';
%logDir = '../../data/mhfastslam/';
logDir = '/datasets/RBPHDSLAM/Pd0.99_c0.0001/rbphdslam/';
logDir = '/datasets/RBPHDSLAM/Pd0.99_c0.5000/rbphdslam/';
logDir = '/datasets/RBPHDSLAM/Pd0.50_c0.0001/rbphdslam/';
logDir = '/datasets/RBPHDSLAM/Pd0.50_c0.5000/rbphdslam/';
%logDir = '/datasets/RBPHDSLAM/Pd0.99_c0.0001/fastslam/';
%logDir = '/datasets/RBPHDSLAM/Pd0.99_c0.5000/fastslam/';
%logDir = '/datasets/RBPHDSLAM/Pd0.50_c0.0001/fastslam/';
%logDir = '/datasets/RBPHDSLAM/Pd0.50_c0.5000/fastslam/';
%logDir = '/datasets/RBPHDSLAM/Pd0.99_c0.0001/mhfastslam/';
%logDir = '/datasets/RBPHDSLAM/Pd0.99_c0.5000/mhfastslam/';
%logDir = '/datasets/RBPHDSLAM/Pd0.50_c0.0001/mhfastslam/';
%logDir = '/datasets/RBPHDSLAM/Pd0.50_c0.5000/mhfastslam/';
showLastFrameOnly = 1;
showGroundTruthRobot = 0;
showEstimates = 1;

disp('Opening ground-truth pose file');
fid_gtPose = fopen(strcat(logDir, 'gtPose.dat'));
disp('Opening ground-truth landmark file');
fid_gtLmk = fopen(strcat(logDir, 'gtLandmark.dat'));
disp('Opening odometry file');
fid_odo = fopen(strcat(logDir, 'odometry.dat'));
disp('Opening dead-reckoning file');
fid_dr = fopen(strcat(logDir, 'deadReckoning.dat'));
disp('Opening measurement file');
fid_meas = fopen(strcat(logDir, 'measurement.dat'));
disp('Opening particle pose file');
fid_pEst = fopen(strcat(logDir, 'particlePose.dat'));
disp('Opening landmark estimate file');
fid_lmEst = fopen(strcat(logDir, 'landmarkEst.dat'));

gt_pose = fscanf(fid_gtPose, '%f %f %f %f\n', [4, inf]);
gt_lmk = fscanf(fid_gtLmk, '%f %f %f\n', [3, inf]);
dr_pose = fscanf(fid_dr, '%f %f %f %f\n', [4, inf]);
z_k = textscan(fid_meas, '%f %f %f', 1);  

kMax = textscan(fid_pEst, 'Timesteps: %f');
%nParticles = textscan(fid_pEst, 'nParticles: %d\n');
kMax = textscan(fid_lmEst, 'Timesteps: %f');
nParticles = textscan(fid_lmEst, 'nParticles: %d\n');
%odom = fscanf(fid_odo, '%f %f %f %f\n', [4, inf]);
%meas = fscanf(fid_meas, '%f %f %f\n', [3, inf]);

% Initiate plot
hfig = figure('Renderer','OpenGL');
set(gcf, 'Color', 'w');
set(hfig,'Position',[0,0,1024,1024*3/4]);
title('Groundturth robot trajectory and landmark positions');
if(showGroundTruthRobot)
    plot(gt_pose(2,:), gt_pose(3,:), 'r--');
    hold on
end
plot(gt_lmk(1,:), gt_lmk(2,:), 'k.','MarkerSize',7);
hold on
axis equal
grid on
set(gca, 'XLim', get(gca, 'XLim') + [-1, 1] );
set(gca, 'YLim', get(gca, 'YLim') + [-1, 1] );
tmp = get(gca, 'XLim');
text_x = (tmp(2) - tmp(1)) * 9/10 + tmp(1);
tmp = get(gca, 'YLim');
text_y = (tmp(2) - tmp(1)) * 9.5/10 + tmp(1);

est_pose = zeros(3, kMax{1});

h_particlePos = 0;
h_time = text(text_x, text_y, sprintf('%d',0));
if(showGroundTruthRobot)
    h_robotPos = plot(gt_pose(2,1), gt_pose(3,1), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    h_robotHdg = line([gt_pose(2,1) gt_pose(2,1)+0.5*cos(gt_pose(4,1))], [gt_pose(3,1) gt_pose(3,1)+0.5*sin(gt_pose(4,1))], 'Color', 'k');
end
%h_drPos = plot(dr_pose(2,1), dr_pose(3,1), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
%h_drHdg = line([dr_pose(2,1) dr_pose(2,1)+0.5*cos(dr_pose(4,1))], [dr_pose(3,1) dr_pose(3,1)+0.5*sin(dr_pose(4,1))], 'Color', 'k');

for k = 1:kMax{1} % actual time index starts at 0
    
    kSim = textscan(fid_pEst, 'k = %d\n');
    
    % Refresh plot
    if(showLastFrameOnly == 0 || k == 1)
            
        if(showGroundTruthRobot)
            delete( h_robotPos )
            delete( h_robotHdg )
        end
        %delete( h_drPos )
        %delete( h_drHdg )
        if h_particlePos ~= 0
            delete( h_particlePos );
        end
        delete(findobj('Color',[0.4, 0.4, 1]))
        delete(findobj('Color','g'))
        delete(findobj('LineWidth',1.1))
        delete( h_time )
    end
    
    %Plot groundtruth pose
    
    if(showLastFrameOnly == 0 || (showLastFrameOnly == 1 && k == kMax{1}) )
        t = gt_pose(1,k);
        x = gt_pose(2,k);
        y = gt_pose(3,k);
        z = gt_pose(4,k);
        if(showGroundTruthRobot)    
            h_robotPos = plot(x, y, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
            h_robotHdg = line([x x+0.5*cos(z)], [y y+0.5*sin(z)], 'Color', 'k');
        end
    end
    
    %Plot measurements
    while(abs(z_k{1} - (k-1)*0.1) < 1e-10)
        if(showLastFrameOnly == 0 || (showLastFrameOnly == 1 && k == kMax{1}) )
            r = z_k{2}; 
            b = z_k{3};
            mx = x + r*cos(b + z);
            my = y + r*sin(b + z);
            line([x mx], [y my], 'Color', [0.4, 0.4, 1]);
        end
        z_k = textscan(fid_meas, '%f %f %f', 1);
    end
    
    %Plot dead reckoning
    %dr_x = dr_pose(2,k);
    %dr_y = dr_pose(3,k);
    %dr_z = dr_pose(4,k);
    %h_drPos = plot(dr_x, dr_y, 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
    %h_drHdg = line([dr_x dr_x+0.5*cos(dr_z)], [dr_y dr_y+0.5*sin(dr_z)], 'Color', 'g');
    
    %Plot particle pose
    nParticles = textscan(fid_pEst, 'nParticles = %d\n');
    x_i_k = textscan(fid_pEst, '%f %f %f %f\n');

    x_i = zeros(4, nParticles{1});
    highest_weight_i = 0;
    highest_weight = 0;
    for i = 1:nParticles{1}
        x_i(1, i) = x_i_k{1}(i); %x
        x_i(2, i) = x_i_k{2}(i); %y
        x_i(3, i) = x_i_k{3}(i); %theta
        x_i(4, i) = x_i_k{4}(i); %w
        if x_i(4, i) > highest_weight
            highest_weight = x_i(4, i);
            highest_weight_i = i;
        end      
    end
    
    est_pose(1,k) = x_i(1, highest_weight_i);
    est_pose(2,k) = x_i(2, highest_weight_i);
    est_pose(3,k) = x_i(3, highest_weight_i);
        
    if(showEstimates)
        if(showLastFrameOnly == 0 || (showLastFrameOnly == 1 && k == kMax{1}) )
            h_particlePos = plot(x_i(1, :), x_i(2, :), 'm.');
        end
    end
    
    % Plot map for the highest weight particle
    for i = 1:nParticles{1}
        
        header = textscan(fid_lmEst, 'Timestep: %d   Particle: %d   Map Size: %d', 1);
        map_size = header{3};
        data = textscan(fid_lmEst, '%f %f %f %f %f %f %f\n');
    
        if(showEstimates)
        if(showLastFrameOnly == 0 || (showLastFrameOnly == 1 && k == kMax{1}) )

            if i == highest_weight_i
                %display(map_size);
                for m = 1:map_size

                    u = [data{1}(m); data{2}(m)];
                    cov = [data{3}(m) data{4}(m); data{5}(m) data{6}(m)];
                    w = data{7}(m);

                    [evec, eval] = eig(cov);
                    axes_length = 3 * sqrt(diag(eval));
                    if(axes_length(2) > axes_length(1))
                        angle = atan2(evec(2,1), evec(1,1));
                    else
                        angle = atan2(evec(2,2), evec(1,2));
                    end

                    h_ellipse = ellipse(axes_length(1), axes_length(2), angle, u(1), u(2));
                    set(h_ellipse,'color',[max(0.8-w,0),max(0.8-w,0),1]); 
                    set(h_ellipse,'MarkerEdgeColor',[max(0.8-w,0),max(0.8-w,0),1]); 
                    set(h_ellipse,'LineWidth',1.1);

                end
            end
            
        end
        end
        
    end
    
    if(showLastFrameOnly == 0 || (showLastFrameOnly == 1 && k == kMax{1}) )
        h_time = text(text_x, text_y, sprintf('%d',k)); 
    end
    
    if(k == 1)
        mkdir( logDir, 'animation');    
    end
    %export_fig( sprintf(strcat( logDir, '/animation/%06d.png'),k));

    if(showLastFrameOnly == 0 || (showLastFrameOnly == 1 && k == kMax{1}) ) 
        pause(0.002)
    end

end

delete( h_time )
delete(findobj('Color',[0.4, 0.4, 1]))

if(showEstimates)
    plot(est_pose(1,:), est_pose(2,:), 'm--');
end

export_fig( sprintf(strcat( logDir, 'estimate.pdf')));

fclose(fid_gtPose);
fclose(fid_gtLmk);
fclose(fid_odo);
fclose(fid_dr);
fclose(fid_meas);