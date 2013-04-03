% Plotting for PHD Filter Simulation

figure;
title('Vehicle Position');
hold on
plot(1:k_max, p_k_groundtruth(1,:), 'b-');
grid on

for i = 1:n_particles
    plot(k_sim_start:k_sim_end, p_k__i(1, k_sim_start:k_sim_end, i), 'g-' );
    hold on
end
[max_weight, max_index] = max(particle_weight(:,k_max));
plot(k_sim_start:k_sim_end, p_k__i(1, k_sim_start:k_sim_end, max_index), 'r-', 'LineWidth', 2);

figure;
hold on
title('Vehicle Position Estimate Error');
for i = 1:n_particles
    error = p_k__i(1, k_sim_start:k_sim_end, i) - p_k_groundtruth(1, k_sim_start:k_sim_end);
    plot(k_sim_start:k_sim_end, error, 'g-' );
    hold on
end
error = p_k__i(1, k_sim_start:k_sim_end, max_index) - p_k_groundtruth(1, k_sim_start:k_sim_end);
plot(k_sim_start:k_sim_end, error, 'r-', 'LineWidth', 2);
grid on;

x_min = floor( min(p_k__i(1,:)) ) - 5;
x_max = ceil( max(p_k__i(1,:)) ) + 5;
x = x_min : 0.05 : x_max;
y = zeros( 1, length(x) );
figure;
grid on
hold on
title('Feature Hypothesis Density');
for m = 1:M_size(max_index)
    u = M{max_index,1}(1,m);
    cov = M{max_index,2}(1,m);
    w = M{max_index,3}(m);
    y = y + w*pdf('normal', x, u, sqrt(cov));
    plot(u, w, 'k.', 'markersize', 5);
end
plot(x, y, 'r-');
plot(map(1,:), 0.1, 'b.', 'markersize', 5);

figure;
title('timing analysis')
hold on
plot(1:k_max, time_birth, 'r');
plot(1:k_max, time_propogation, 'g');
plot(1:k_max, time_update, 'b');
plot(1:k_max, time_merging, 'k');
plot(1:k_max, time_weighting, 'c');
plot(1:k_max, time_resampling, 'm.');
grid on

figure;
title('particle weights')
hold on
for i = 1:n_particles
    plot(1:k_max, particle_weight(i,1:k_max), 'b-');
end
grid on;
plot(1:k_max, particle_weight(max_index ,1:k_max), 'r-', 'LineWidth', 2);

