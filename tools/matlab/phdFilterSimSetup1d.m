% PHD Filter Simulation setup
% Keith Leung 2013

clear all;
close all;
rand('seed', 6);
randn('seed', 6);

%% Simulation settings
k_max = 1000;
n_features = 10;
y_rangeLim = 5;

force_particle_on_true_trajectory = 0;
no_measurement_noise = 0;

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
noise_motion = [0.03, 1e-15, 1e-15, 1e-15, 1e-15, 1e-15]; 
noise_obs = [0.1, 1e-10, 1e-10];
[p_k, c_k, d_k_noiseFree, D_vec_k_noiseFree] = generate1dTrajectory(k_max, 0.1, [1, 0.1]);
[d_k, D_vec_k] = generate1dOdometry(d_k_noiseFree, D_vec_k_noiseFree, noise_motion);
map = generate1dFeatures(p_k, c_k, n_features, y_rangeLim);
noise_obs_ = noise_obs;
if no_measurement_noise == 1
    noise_obs(1) = noise_obs(1) * 1e-10;
end
[y, Y_vec] = generate1dMeasurements(p_k, c_k, map, y_rangeLim, noise_obs);
noise_obs = noise_obs_;
if no_measurement_noise == 1
    Y_vec(1,:) = Y_vec(1,:) / 1e-10;
end

p_k_groundtruth = p_k;
c_k_groundtruth = c_k;

%% Estimator Settings
%noise_motion = [0.05, 0.05, 0.05, 0.01, 0.01, 0.01];
n_particles = 50;
map_size_limit = 1000; % number of Gaussians in each particle's map
noise_motion = noise_motion * 1; 
noise_obs = noise_obs * 3 ; % inflate noise for Kalman filter
resample_interval = 0;
effective_particle_threshold = n_particles / 2;

% PHD Filter Settings
sensor_limit_upper_buffer = 0.25;
sensor_limit_lower_buffer = 0.25;
particle_weighting_strategy = 1; % 0 for empty strategy, 1 for single feature proper, 2 for single feature without exp term
birth_Gaussian_likelihood_threshold = 2;
merging_mahalanoblis_distance_threshold = 0.25;
feature_merging_covarance_inflation_factor = 2;

% Simulator Seed
rand('seed', 6);
randn('seed', 6);

phdFilterSim1d;
phdFilterSimPlot1d;