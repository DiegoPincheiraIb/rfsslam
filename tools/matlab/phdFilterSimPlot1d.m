% Plotting for PHD Filter Simulation
addpath '~/src/Matlab/export_fig' % export_fig can be downloaded from Matlab Central

figure;
set(gcf, 'Color', 'w');

%title('Vehicle Trajectory');
xlabel('Time [s]');
ylabel('Position [m]');
hold on
% for i = 1:n_particles
%     plot(k_sim_start:k_sim_end, p_k__i(1, k_sim_start:k_sim_end, i), 'b-' );
%     hold on
% end
plot(k_sim_start:k_sim_end, p_k_max_weight(k_sim_start:k_sim_end), 'm-', 'LineWidth', 1);
plot(k_sim_start:k_sim_end, p_k_weighted(k_sim_start:k_sim_end), 'r-', 'LineWidth', 1);
plot(1:k_max, p_k_groundtruth(1,:), 'k-');
plot(1:k_max, p_k_dr(1,:), 'g-');
grid on
%export_fig(sprintf('results/vehicle_trajectory_p%d_ws%d.pdf', n_particles, particle_weighting_feautre_set_max_size ));

%%
figure;
set(gcf, 'Color', 'w');
hold on
xlabel('Time [s]');
ylabel('Position error [m]');
%title('Vehicle Position Estimate Error');
% for i = 1:n_particles
%     error = p_k__i(1, k_sim_start:k_sim_end, i) - p_k_groundtruth(1, k_sim_start:k_sim_end);
%     plot(k_sim_start:k_sim_end, error, 'b-' );
%     hold on
% end
error = p_k_max_weight(k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
plot(k_sim_start:k_sim_end, error, 'm-', 'LineWidth', 1);
error = p_k_weighted(k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
plot(k_sim_start:k_sim_end, error, 'r-', 'LineWidth', 1);
error = p_k_dr(1, k_sim_start:k_sim_end) - p_k_groundtruth(1, k_sim_start:k_sim_end);
plot(k_sim_start:k_sim_end, error, 'g-', 'LineWidth', 1);
grid on;
%export_fig(sprintf('results/vehicle_position_error_p%d_ws%d.pdf', n_particles, particle_weighting_feautre_set_max_size ));


%%

figure;
set(gcf, 'Color', 'w');
grid on
hold on
xlabel('x [m]');
ylabel('Map Intensity');
%title('Feature Hypothesis Density');
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
        %plot(u, w, 'k.', 'markersize', 5);
    end
end
plot(x, y, 'r-');
for m = 1:length(map(1,:))
    line( [map(1,m), map(1,m)], [-1 0], 'Color', 'b', 'LineWidth', 2);
end
ylim([-2 15]);
%export_fig(sprintf('results/map_intensity_p%d_ws%d.pdf', n_particles, particle_weighting_feautre_set_max_size ));

%%
figure;
set(gcf, 'Color', 'w');
grid on;
hold on;
xlabel('Time [s]');
ylabel('Feature Count');
%title('Map Dimensionality Estimate');
plot(k_sim_start:k_sim_end, nFeaturesEstimateAllParticles(k_sim_start:k_sim_end), 'r-');
%plot(k_sim_start:k_sim_end, nFeaturesEstimateMaxWeightParticle(k_sim_start:k_sim_end), 'm-');
plot(k_sim_start:k_sim_end, nFeaturesObserved(k_sim_start:k_sim_end), 'k-', 'LineWidth', 2);
ylim([0 15]);
%export_fig(sprintf('results/map_count_p%d_ws%d.pdf', n_particles, particle_weighting_feautre_set_max_size ));


%%
figure;
set(gcf, 'Color', 'w');
grid on;
hold on;
xlabel('Time [s]');
ylabel('Map Error [m]');
%title('Map Estimate Error');
plot(k_sim_start:k_sim_end, map_estimate_error_all_particles(k_sim_start:k_sim_end), 'r-', 'LineWidth', 1);
%plot(k_sim_start:k_sim_end, map_estimate_error_highest_weight_particle(k_sim_start:k_sim_end), 'm-', 'LineWidth', 2);
ylim([0 5])
%export_fig(sprintf('results/map_error_p%d_ws%d.pdf', n_particles, particle_weighting_feautre_set_max_size ));


% figure;
% title('timing analysis')
% hold on
% plot(1:k_max, time_birth, 'r');
% plot(1:k_max, time_propogation, 'g');
% plot(1:k_max, time_update, 'b');
% plot(1:k_max, time_merging, 'k');
% plot(1:k_max, time_weighting, 'c');
% plot(1:k_max, time_resampling, 'm.');
% grid on

% figure;
% title('particle weights')
% hold on
% for i = 1:n_particles
%     plot(1:k_max, particle_weight(i,1:k_max), 'b-');
% end
% grid on;
% plot(1:k_max, particle_weight(max_index ,1:k_max), 'r-', 'LineWidth', 2);

