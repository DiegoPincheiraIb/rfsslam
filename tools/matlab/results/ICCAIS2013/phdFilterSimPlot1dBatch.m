clear all
close all

load multiFeature_Nc025.mat

%% Robot trajectory
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [220 540 560 420]);
xlabel('Time [s]');
ylabel('x-position [m]');
hold on
grid on
plot(1:k_max, p_k_groundtruth(1,:), 'k-', 'LineWidth', 3);
plot(1:k_max, p_k_dr(1,:), 'g-');
h_robotTrajectory = gca;
h_robotTrajectotyFig = gcf;

%% Position Estimate Error
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [220 40 560 420]);
xlabel('Time [s]');
ylabel('x-position error error [m]');
hold on
grid on
error = p_k_dr(1, k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
plot(k_sim_start:k_sim_end, error, 'g-', 'LineWidth', 1);
h_estimateError = gca;
h_estimateErrorFig = gcf;

%% Estimated Feature Count
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [780 540 560 420]);
grid on;
hold on;
xlabel('Time [s]');
ylabel('Estimated Feature Count');
plot(k_sim_start:k_sim_end, nFeaturesObserved(k_sim_start:k_sim_end), 'k-', 'LineWidth', 3);
ylim([0 25]);
h_estimatedFeatureCount = gca;
h_estimatedFeatureCountFig = gcf;

%% Map Error
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [780 40 560 420]);
grid on;
hold on;
xlabel('Time [s]');
ylabel('Map Estimate Error [m]');
ylim([0 5])
h_mapEstimateError = gca;
h_mapEstimateErrorFig = gcf;

%% Plot estimator results

%% Plot multi-feature strategy estimation results
plot(h_robotTrajectory, k_sim_start:k_sim_end, p_k_max_weight(k_sim_start:k_sim_end), 'r-', 'LineWidth', 1);
%plot(h_robotTrajectory, k_sim_start:k_sim_end, p_k_weighted(k_sim_start:k_sim_end), 'r-', 'LineWidth', 1);
error = p_k_max_weight(k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
plot(h_estimateError, k_sim_start:k_sim_end, error, 'r-', 'LineWidth', 1);
%error = p_k_weighted(k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
%plot(h_estimatError, k_sim_start:k_sim_end, error, 'r-', 'LineWidth', 1);
plot(h_estimatedFeatureCount, k_sim_start:10:k_sim_end, nFeaturesEstimateAllParticles(k_sim_start:10:k_sim_end), 'r-', 'LineWidth', 2);
%plot(h_estimatedFeatureCount, k_sim_start:k_sim_end, nFeaturesEstimateMaxWeightParticle(k_sim_start:k_sim_end), 'r-');
plot(h_mapEstimateError, k_sim_start:20:k_sim_end, map_estimate_error_all_particles(k_sim_start:20:k_sim_end), 'r-', 'LineWidth', 2);
%plot(h_mapEstimateError, k_sim_start:5:k_sim_end, map_estimate_error_highest_weight_particle(k_sim_start:5:k_sim_end), 'r-', 'LineWidth', 1);

% Map intensity
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [1350 10 560 320])
xlabel('x-position [m]');
ylabel('Map Intensity');
grid on
hold on
total_particle_weight = sum(particle_weight(:,k_max));
x_min = floor( min(p_k__i(1,:)) ) - 5;
x_max = ceil( max(p_k__i(1,:)) ) + 5;
x = x_min : 0.025 : x_max;
y = zeros( 1, length(x) );
for i = 1:n_particles
    for m = 1:M_size(i)
        u = M{i,1}(1,m);
        cov = M{i,2}(1,m);
        w = M{i,3}(m);
        y = y + w*pdf('normal', x, u, sqrt(cov))*particle_weight(i,k_max)/total_particle_weight;
    end
end
plot(x, y, 'r-');
for m = 1:length(map(1,:))
    line( [map(1,m), map(1,m)], [-1 0], 'Color', 'b', 'LineWidth', 2);
end
xlim([-5 30]);
ylim([-2 15]);
h_mapMultiFeature = gca;
h_mapMultiFeatureFig = gcf;
legend('v^+(x)', 'Landmark Positions')

%% Plot empty-set strategy estimation results
load singleFeature_Nc025.mat
plot(h_robotTrajectory, k_sim_start:k_sim_end, p_k_max_weight(k_sim_start:k_sim_end), 'b-', 'LineWidth', 1);
%plot(h_robotTrajectory, k_sim_start:k_sim_end, p_k_weighted(k_sim_start:k_sim_end), 'b-', 'LineWidth', 1);
error = p_k_max_weight(k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
plot(h_estimateError, k_sim_start:k_sim_end, error, 'b-', 'LineWidth', 1);
%error = p_k_weighted(k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
%plot(h_estimatError, k_sim_start:k_sim_end, error, 'b-', 'LineWidth', 1);
plot(h_estimatedFeatureCount, k_sim_start:10:k_sim_end, nFeaturesEstimateAllParticles(k_sim_start:10:k_sim_end), 'b-');
%plot(h_estimatedFeatureCount, k_sim_start:k_sim_end, nFeaturesEstimateMaxWeightParticle(k_sim_start:k_sim_end), 'b-');
plot(h_mapEstimateError, k_sim_start:20:k_sim_end, map_estimate_error_all_particles(k_sim_start:20:k_sim_end), 'b-', 'LineWidth', 1);
%plot(h_mapEstimateError, k_sim_start:5:k_sim_end, map_estimate_error_highest_weight_particle(k_sim_start:5:k_sim_end), 'b-', 'LineWidth', 1);

% Map intensity
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [1350 330 560 320])
xlabel('x-position [m]');
ylabel('Map Intensity');
grid on
hold on
total_particle_weight = sum(particle_weight(:,k_max));
x_min = floor( min(p_k__i(1,:)) ) - 5;
x_max = ceil( max(p_k__i(1,:)) ) + 5;
x = x_min : 0.025 : x_max;
y = zeros( 1, length(x) );
for i = 1:n_particles
    for m = 1:M_size(i)
        u = M{i,1}(1,m);
        cov = M{i,2}(1,m);
        w = M{i,3}(m);
        y = y + w*pdf('normal', x, u, sqrt(cov))*particle_weight(i,k_max)/total_particle_weight;
    end
end
plot(x, y, 'r-');
for m = 1:length(map(1,:))
    line( [map(1,m), map(1,m)], [-1 0], 'Color', 'b', 'LineWidth', 2);
end
xlim([-5 30]);
ylim([-2 15]);
h_mapSingleFeature = gca;
h_mapSingleFeatureFig = gcf;
legend('v^+(x)', 'Landmark Positions')

%% Plot empty-set strategy estimation results
load emptySet_Nc025.mat
plot(h_robotTrajectory, k_sim_start:k_sim_end, p_k_max_weight(k_sim_start:k_sim_end), 'm-', 'LineWidth', 1);
%plot(h_robotTrajectory, k_sim_start:k_sim_end, p_k_weighted(k_sim_start:k_sim_end), 'm-', 'LineWidth', 1);
error = p_k_max_weight(k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
plot(h_estimateError, k_sim_start:k_sim_end, error, 'm-', 'LineWidth', 1);
%error = p_k_weighted(k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
%plot(h_estimatError, k_sim_start:k_sim_end, error, 'm-', 'LineWidth', 1);
plot(h_estimatedFeatureCount, k_sim_start:10:k_sim_end, nFeaturesEstimateAllParticles(k_sim_start:10:k_sim_end), 'm-');
%plot(h_estimatedFeatureCount, k_sim_start:k_sim_end, nFeaturesEstimateMaxWeightParticle(k_sim_start:k_sim_end), 'm-');
plot(h_mapEstimateError, k_sim_start:20:k_sim_end, map_estimate_error_all_particles(k_sim_start:20:k_sim_end), 'm-', 'LineWidth', 1);
%plot(h_mapEstimateError, k_sim_start:5:k_sim_end, map_estimate_error_highest_weight_particle(k_sim_start:5:k_sim_end), 'm-', 'LineWidth', 1);

% Map intensity
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [1350 650 560 320])
xlabel('x-position [m]');
ylabel('Map Intensity');
grid on
hold on
total_particle_weight = sum(particle_weight(:,k_max));
x_min = floor( min(p_k__i(1,:)) ) - 5;
x_max = ceil( max(p_k__i(1,:)) ) + 5;
x = x_min : 0.025 : x_max;
y = zeros( 1, length(x) );
for i = 1:n_particles
    for m = 1:M_size(i)
        u = M{i,1}(1,m);
        cov = M{i,2}(1,m);
        w = M{i,3}(m);
        y = y + w*pdf('normal', x, u, sqrt(cov))*particle_weight(i,k_max)/total_particle_weight;
    end
end
plot(x, y, 'r-');
for m = 1:length(map(1,:))
    line( [map(1,m), map(1,m)], [-1 0], 'Color', 'b', 'LineWidth', 2);
end
xlim([-5 30]);
ylim([-2 15]);
h_mapEmptySet = gca;
h_mapEmptySetFig = gcf;
legend('v^+(x)', 'Landmark Positions')

%% Plot MH-FastSLAM results

load FastSLAM_Nc025.mat
plot(h_robotTrajectory, k_sim_start:k_sim_end, p_k_max_weight(k_sim_start:k_sim_end), 'k--', 'LineWidth', 1);
%plot(h_robotTrajectory, k_sim_start:k_sim_end, p_k_weighted(k_sim_start:k_sim_end), 'm-', 'LineWidth', 1);
error = p_k_max_weight(k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
plot(h_estimateError, k_sim_start:k_sim_end, error, 'k--', 'LineWidth', 1);
%error = p_k_weighted(k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
%plot(h_estimatError, k_sim_start:k_sim_end, error, 'm-', 'LineWidth', 1);
plot(h_estimatedFeatureCount, k_sim_start:10:k_sim_end, nFeaturesEstimateAllParticles(k_sim_start:10:k_sim_end), 'k--');
%plot(h_estimatedFeatureCount, k_sim_start:k_sim_end, nFeaturesEstimateMaxWeightParticle(k_sim_start:k_sim_end), 'm-');
plot(h_mapEstimateError, k_sim_start:20:k_sim_end, map_estimate_error_all_particles(k_sim_start:20:k_sim_end), 'k--', 'LineWidth', 1);
%plot(h_mapEstimateError, k_sim_start:5:k_sim_end, map_estimate_error_highest_weight_particle(k_sim_start:5:k_sim_end), 'm-', 'LineWidth', 1);

% Map intensity
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [1350 650 560 320])
xlabel('x-position [m]');
ylabel('Map Intensity');
grid on
hold on
total_particle_weight = sum(particle_weight(1:n_particles,k_max));
x_min = floor( min(p_k__i(1,:)) ) - 5;
x_max = ceil( max(p_k__i(1,:)) ) + 5;
x = x_min : 0.025 : x_max;
y = zeros( 1, length(x) );
for i = 1:n_particles
    for m = 1:M_size(i)
        u = M{i,1}(1,m);
        cov = M{i,2}(1,m);
        log_w = M{i,3}(m);
        w = 1 - 1/(1 + exp(log_w));
        y = y + w*pdf('normal', x, u, sqrt(cov))*particle_weight(i,k_max)/total_particle_weight;
    end
end
plot(x, y, 'r-');
for m = 1:length(map(1,:))
    line( [map(1,m), map(1,m)], [-1 0], 'Color', 'b', 'LineWidth', 2);
end
xlim([-5 30]);
ylim([-2 15]);
h_mapFastSLAM = gca;
h_mapFastSLAMFig = gcf;
legend('v^+(x)', 'Landmark Positions')


%%
legend(h_robotTrajectory,'Groundtruth','Dead reckoning', 'Multi-feature strategy', 'Single-feature strategy', 'Empty-set strategy', 'MH-FastSLAM', 'Location', 'NorthWest')
legend(h_estimateError,'Dead reckoning', 'Multi-feature strategy', 'Single-feature strategy', 'Empty-set strategy', 'MH-FastSLAM', 'Location', 'NorthWest')
legend(h_estimatedFeatureCount,'Groundtruth', 'Multi-feature strategy', 'Single-feature strategy', 'Empty-set strategy', 'MH-FastSLAM', 'Location', 'SouthEast')
legend(h_mapEstimateError,'Multi-feature strategy', 'Single-feature strategy', 'Empty-set strategy', 'MH-FastSLAM', 'Location', 'NorthWest')

%%

return;

addpath '~/src/Matlab/export_fig'
export_fig('robotTrajectoryEstimate.pdf', h_robotTrajectory);
export_fig('robotTrajectoryEstimateError.pdf', h_estimateError);
export_fig('featureCountEstimate.pdf', h_estimatedFeatureCount);
export_fig('mapError.pdf', h_mapEstimateError);
export_fig('mapMultiFeature.pdf', h_mapMultiFeature);
export_fig('mapSingleFeature.pdf', h_mapSingleFeature);
export_fig('mapEmptySet.pdf', h_mapEmptySet);
export_fig('mapFastSLAM.pdf', h_mapFastSLAM);

saveas(h_robotTrajectotyFig, 'robotTrajectoryEstimate.fig', 'fig');
saveas(h_estimateErrorFig, 'robotTrajectoryEstimateError.fig', 'fig');
saveas(h_estimatedFeatureCountFig, 'featureCountEstimate.fig', 'fig');
saveas(h_mapEstimateErrorFig, 'mapError.fig', 'fig');
saveas(h_mapMultiFeatureFig, 'mapMultiFeature.fig', 'fig');
saveas(h_mapSingleFeatureFig, 'mapSingleFeature.fig', 'fig');
saveas(h_mapEmptySetFig, 'mapEmptySet.fig', 'fig');
saveas(h_mapFastSLAMFig, 'mapFastSLAM.fig', 'fig');
