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

disp('Reading particle pose file');
fid = fopen('../../data/particlePose.dat');
kMax = textscan(fid, 'Timesteps: %d');
nParticles = textscan(fid, 'nParticles: %d\n');
x_i = zeros(3, nParticles{1}, kMax{1});
for k = 1:kMax{1} 
    kSim = textscan(fid, 'k = %d\n');
    x_i_temp = textscan(fid, '%f %f %f\n');
    for i = 1:nParticles{1}
        x_i(1, i, k) = x_i_temp{1}(i);
        x_i(2, i, k) = x_i_temp{2}(i);
        x_i(3, i, k) = x_i_temp{3}(i);
    end
end

disp('Reading landmark estimate file');
fid = fopen('../../data/landmarkEst.dat');
kMax = textscan(fid, 'Timesteps: %d');
nParticles = textscan(fid, 'nParticles: %d\n');
landmarkEst = cell(nParticles{1}, kMax{1} );
for k = 1:kMax{1} - 1
    for i = 1:nParticles{1}
        header = textscan(fid, 'Timestep: %d   Particle: %d   Map Size: %d', 1);
        map_size = header{3};
        data = textscan(fid, '%f %f %f %f %f %f %f\n');
        landmarkEst{i, k + 1} = cell(map_size, 3);
        for m = 1:map_size
            u = [data{1}(m); data{2}(m)];
            cov = [data{3}(m) data{4}(m); data{5}(m) data{6}(m)];
            weight = data{7}(m);
            landmarkEst{i, k + 1}{m, 1} = u;
            landmarkEst{i, k + 1}{m, 2} = cov;
            landmarkEst{i, k + 1}{m, 3} = weight;
        end
    end
end

playCPP2dSimData;
   
