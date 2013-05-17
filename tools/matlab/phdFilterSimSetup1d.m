% PHD Filter Simulation setup
% Keith Leung 2013

clear all;
close all;

%% Simulation settings
k_max = 3000;
n_features = 10;
y_rangeLim = 5;
force_particle_on_true_trajectory = 0;
no_measurement_noise = 0;
simulator_seed = 6;

%% Simulation Data Settings

noise_motion = [0.1, 1e-15, 1e-15, 1e-15, 1e-15, 1e-15]; 
noise_obs = [0.1, 1e-10, 1e-10];
P_detection_static = 0.95;
N_c = 1; % number of expected clutter measurements per timestep
clutter_intensity = N_c / (y_rangeLim * 2); % clutter intensity
data_generation_seed = 18; %8-8, 68-8, 18-10, 18-np=10-k=3000

%% Estimator Settings
n_particles = 25;
map_size_limit = 1000; % number of Gaussians in each particle's map
noise_motion_inflation_factor = 1; 
noise_obs_inflation_factor = 2 ; % inflate noise for Kalman filter

resample_interval = 0;
effective_particle_threshold = n_particles / 4;

sensor_limit_upper_buffer = 0.20;
sensor_limit_lower_buffer = 0.20;
particle_weighting_feautre_set_max_size = 0;
particle_weighting_random_map = 0;
birth_Gaussian_likelihood_threshold = 2;
merging_mahalanoblis_distance_threshold = 0.1;
merging_weight_threshold = 2;
feature_merging_covarance_inflation_factor = 5;

%% Evaluation Settings
cutoff = 3; 

%% Visualization Settings
visualize = 0;

%% NO MORE CONFIGURATIONS FROM HERE ON 

% Generate Simulation Data

rand('seed', data_generation_seed); 
randn('seed', data_generation_seed); 

[p_k, c_k, d_k_noiseFree, D_vec_k_noiseFree] = generate1dTrajectory(k_max, 0.1, [1, 0.1]);
[d_k, D_vec_k] = generate1dOdometry(d_k_noiseFree, D_vec_k_noiseFree, noise_motion);
map = generate1dFeatures(p_k, c_k, n_features, y_rangeLim);
noise_obs_ = noise_obs;
if no_measurement_noise == 1
    noise_obs(1) = noise_obs(1) * 1e-10;
end
[y, Y_vec] = generate1dMeasurements(p_k, c_k, map, y_rangeLim, noise_obs, P_detection_static, N_c);
P_detection_static = sum(y(5,:) ~= -1)/length(y(5,:));
noise_obs = noise_obs_;
if no_measurement_noise == 1
    Y_vec(1,:) = Y_vec(1,:) / 1e-10;
end

p_k_groundtruth = p_k;
c_k_groundtruth = c_k;

% Apply simulator settings

noise_motion = noise_motion * noise_motion_inflation_factor;
noise_obs = noise_obs * noise_obs_inflation_factor;

rand('seed', simulator_seed);
randn('seed', simulator_seed);

phdFilterSim1d;
phdFilterSimPlot1d;