% Plotting for PHD Filter Simulation

figure;
hold on
plot3(p_k_groundtruth(1,:), p_k_groundtruth(2,:), p_k_groundtruth(3,:), 'b-');
plot3(map(1,:), map(2,:), map(3,:), 'k.');
grid on
axis equal

% plot trajectories of all particles and find the one closest to the
% groundtruth

d_min = norm(p_k_groundtruth(:, k_sim_end) - p_k__i(:, k_sim_end, 1) );
min_dist_index = 1;
min_dist_weight = particle_weight(1, k_max);
for i = 1:n_particles
    plot3(p_k__i(1, k_sim_start:k_sim_end, i), p_k__i(2, k_sim_start:k_sim_end, i), p_k__i(3, k_sim_start:k_sim_end, i), 'g-' );
    hold on
    
    d = norm(p_k_groundtruth(:, k_sim_end) - p_k__i(:, k_sim_end, i) );
    if d < d_min
        d_min = d;
        min_dist_index = i;
        min_dist_weight = particle_weight(i, k_max);
    end
end

[max_weight, max_index] = max(particle_weight(:,k_max));
display(max_weight);
display(min_dist_weight);

plot3(p_k__i(1, k_sim_start:k_sim_end, max_index), p_k__i(2, k_sim_start:k_sim_end, max_index), p_k__i(3, k_sim_start:k_sim_end, max_index), 'r-', 'LineWidth', 2);
for m = 1:M_size(max_index)
    x = M{max_index,1}(:,m);
    if M{max_index,3}(m) > 0.25
        plot3(x(1), x(2), x(3), 'ro'); 
    end
end
plot3(p_k__i(1, k_sim_start:k_sim_end, min_dist_index), p_k__i(2, k_sim_start:k_sim_end, min_dist_index), p_k__i(3, k_sim_start:k_sim_end, min_dist_index), 'm-', 'LineWidth', 2);
for m = 1:M_size(min_dist_index)
    x = M{min_dist_index,1}(:,m);
    if M{min_dist_index,3}(m) > 0.25
        plot3(x(1), x(2), x(3), 'mo'); 
    end
end
grid on
axis equal

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

%figure;
%hold on
%plot(1:k_max, map_size_birth, 'r');
%plot(1:k_max, map_size_update, 'b');
%plot(1:k_max, map_size_merging, 'k');
%grid on

figure;
title('particle weights')
hold on
for i = 1:n_particles
    plot(1:k_max, particle_weight(i,1:k_max), 'b-');
end
grid on;
plot(1:k_max, particle_weight(max_index ,1:k_max), 'r-', 'LineWidth', 2);
plot(1:k_max, particle_weight(min_dist_index ,1:k_max), 'm-', 'LineWidth', 2);

% figure;
% title('measurement likelihood factor')
% hold on
% for i = 1:n_particles
%     plot(1:k_max, measurement_likelihood_factor(i,:), 'b-');
% end
% plot(1:k_max, measurement_likelihood_factor(max_index ,:), 'r-', 'LineWidth', 2);
% grid on
% 
% figure;
% title('similarity factor')
% hold on
% for i = 1:n_particles
%     plot(1:k_max, similarity_factor(i,:), 'b-');
% end
% plot(1:k_max, similarity_factor(max_index ,:), 'r-', 'LineWidth', 2);
% grid on
% 
% figure;
% title('feature count factor')
% hold on
% for i = 1:n_particles
%     plot(1:k_max, feature_count_factor(i,:), 'b-');
% end
% plot(1:k_max, feature_count_factor(max_index ,:), 'r-', 'LineWidth', 2);
% grid on
