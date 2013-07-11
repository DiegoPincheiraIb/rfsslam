% Read the data files generated by c++ program sim2d 
% located in [project home]/data/
% Filenames: gtPose.dat, gtLandmark.dat, odometry.dat, measurement.dat
% Keith Leung 2013

addpath '~/src/Matlab/export_fig'

clear all
close all

disp('Reading ground-truth pose file');
fid = fopen('../../data/gtPose.dat');
gt_pose = fscanf(fid, '%f %f %f %f\n', [4, inf]);
fclose(fid);

disp('Reading ground-truth landmark file');
fid = fopen('../../data/gtLandmark.dat');
gt_lmk = fscanf(fid, '%f %f\n', [2, inf]);
fclose(fid);

disp('Reading odometry file');
fid = fopen('../../data/odometry.dat');
odom = fscanf(fid, '%f %f %f %f\n', [4, inf]);
fclose(fid);

disp('Reading dead-reckoning file');
fid = fopen('../../data/deadReckoning.dat');
dr_pose = fscanf(fid, '%f %f %f %f\n', [4, inf]);
fclose(fid);

disp('Reading measurement file');
fid = fopen('../../data/measurement.dat');
meas = fscanf(fid, '%f %f %f\n', [3, inf]);
fclose(fid);

hfig = figure;
set(gcf, 'Color', 'w');
title('Groundturth robot trajectory and landmark positions');
plot(gt_pose(2,:), gt_pose(3,:), 'r--');
hold on
plot(gt_lmk(1,:), gt_lmk(2,:), 'k.');
axis square
grid on
set(gca, 'XLim', get(gca, 'XLim') + [-1, 1] );
set(gca, 'YLim', get(gca, 'YLim') + [-1, 1] );
   
h_robotPos = plot(gt_pose(2,1), gt_pose(3,1), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
h_robotHdg = line([gt_pose(2,1) gt_pose(2,1)+0.5*cos(gt_pose(4,1))], [gt_pose(3,1) gt_pose(3,1)+0.5*sin(gt_pose(4,1))], 'Color', 'k');
h_drPos = plot(dr_pose(2,1), dr_pose(3,1), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
h_drHdg = line([dr_pose(2,1) dr_pose(2,1)+0.5*cos(dr_pose(4,1))], [dr_pose(3,1) dr_pose(3,1)+0.5*sin(dr_pose(4,1))], 'Color', 'k');
pause(0.005);
meas_idx = 1;
for k = 1 : length(gt_pose)

    delete( h_robotPos )
    delete( h_robotHdg )
    delete( h_drPos )
    delete( h_drHdg )
    delete(findobj('Color','b'))
    
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
    
    
    while( meas(1, meas_idx) == k )
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
    
    export_fig(sprintf('results/anim/%06d.png',k), hfig);
    %pause(0.005)
    
end
    
   
