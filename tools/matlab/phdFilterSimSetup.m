% PHD Filter Simulation setup
% Keith Leung 2013

clear all;
close all;
rand('seed', 3);
randn('seed', 3);

%% Simulation settings
k_max = 1000;
n_features = 5;
y_rangeLim = 5;

%% Generate Simulation Data

% 3d Random Trajectory
% noise_motion = [0.005, 0.005, 0.005, 0.002, 0.002, 0.002]; 
% noise_obs = [0.05, 0.01, 0.01];
%[p_k, c_k, d_k_noiseFree, D_vec_k_noiseFree] = generateTrajectory(k_max, 0.1, [1, 0.1]);
% [d_k, D_vec_k] = generateOdometry(d_k_noiseFree, D_vec_k_noiseFree, noise_motion);
% [y, Y_vec, map] = generateMeasurements(p_k, c_k, n_features, y_rangeLim, noise_obs);

% 3d Circular Trajectory
% noise_motion = [0.005, 0.005, 0.005, 0.002, 0.002, 0.002]; 
% noise_obs = [0.05, 0.01, 0.01];
% [p_k, c_k, d_k_noiseFree, D_vec_k_noiseFree] = generateCircleTrajectory(k_max, 0.1, [1, 0.1]);
% [d_k, D_vec_k] = generateOdometry(d_k_noiseFree, D_vec_k_noiseFree, noise_motion);
% [y, Y_vec, map] = generateMeasurements(p_k, c_k, n_features, y_rangeLim, noise_obs);

% 1d Trajectory
noise_motion = 10*[0.005, 1e-15, 1e-15, 1e-15, 1e-15, 1e-15]; 
noise_obs = [1, 1e-10, 1e-10];
[p_k, c_k, d_k_noiseFree, D_vec_k_noiseFree] = generate1dTrajectory(k_max, 0.1, [1, 0.1]);
[d_k, D_vec_k] = generate1dOdometry(d_k_noiseFree, D_vec_k_noiseFree, noise_motion);
[y, Y_vec, map] = generate1dMeasurements(p_k, c_k, n_features, y_rangeLim, noise_obs);

p_k_groundtruth = p_k;
c_k_groundtruth = c_k;

%% Estimator Settings
%noise_motion = [0.05, 0.05, 0.05, 0.01, 0.01, 0.01];
n_particles = 20;
map_size_limit = 1000; % number of Gaussians in each particle's map
noise_motion = noise_motion * 1; 
noise_obs = noise_obs * 3; % inflate noise for Kalman filter
resample_interval = 0;
particle_weighting_strategy = 1; % 0 for empty strategy, 1 for single feature

% figure
% title('Dead reckoning Monte-Carlo simulation')
% hold on
% plot3(p_k_groundtruth(1,:), p_k_groundtruth(2,:), p_k_groundtruth(3,:), 'r-');
% for n = 1:100
%     [d_k, D_vec_k] = generateOdometry(d_k_noiseFree, D_vec_k_noiseFree, noise_motion);
%     [p_k_dr, C_k_dr] = deadReckoning(d_k, D_vec_k);
%     plot3(p_k_dr(1,:), p_k_dr(2,:), p_k_dr(3,:), 'b-');
% end
% grid on
% axis equal