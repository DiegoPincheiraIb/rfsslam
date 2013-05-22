% MH-FastSLAM Filter 
% Keith Leung 2013

% Vehicle trajectory for each particle
p_k__i = zeros(3, k_max, n_particles); 
C_vec_k__i = zeros(9, k_max, n_particles); 
for i = 1:n_particles
    C_vec_k__i(:, 1, i) = reshape(eye(3), 9, 1); 
end
p_k_dr = zeros(3, k_max);
C_k_dr = zeros(9, k_max);
C_k_dr(:,1) = reshape(eye(3), 9, 1); 
p_k_weighted = zeros(1, k_max);
p_k_max_weight = zeros(1, k_max);

% Per particle map
M = cell(n_particles, 3); % M{:,1} mean, M{:,2} cov reshaped as 9x1 vector, M{:,3} weight 
M_size = zeros(n_particles, 1); % for keeping track of map size
for i = 1:n_particles
    M{i,1} = zeros(3, map_size_limit); % pre-allocate memory for map limit
    M{i,2} = zeros(9, map_size_limit); 
    M{i,3} = zeros(1, map_size_limit);
end
featuresObserved = zeros(n_features, 1);
nFeaturesObserved = zeros(k_max, 1);
nFeaturesEstimate = zeros(n_particles, k_max);
nFeaturesEstimateAllParticles = zeros(k_max, 1);
nFeaturesEstimateMaxWeightParticle = zeros(k_max, 1);
map_estimate_error_all_particles = zeros(k_max, 1);
map_estimate_error_highest_weight_particle = zeros(k_max, 1);

% Particle weighting
particle_weight = ones(n_particles, k_max);
w_before_update = zeros(n_particles, k_max); 
w_after_update = zeros(n_particles, k_max); 
v_before_update = zeros(n_particles, k_max); % For single feature strategy
v_after_update = zeros(n_particles, k_max); % For single feature strategy
measurement_likelihood_factor = zeros(n_particles, k_max);
similarity_factor = zeros(n_particles, k_max);
feature_count_factor = zeros(n_particles, k_max);

% Performance evaluation
time_birth = zeros(1, k_max);
time_propogation = zeros(1, k_max);
time_update = zeros(1, k_max);
time_sorting = zeros(1, k_max);
time_merging = zeros(1, k_max);
time_weighting = zeros(1, k_max);
time_resampling = zeros(1, k_max);
map_size_birth = zeros(1, k_max);
map_size_update = zeros(1, k_max);
map_size_merging = zeros(1, k_max);

%% Some book-keeping
birth_Gaussian_mDist_min = []; % This is used to reduce number of birth Gaussians

k_first_measurement = y(1,1);
idx_prev_obs_start= 1;
idx_current_obs_start= 1;
% find the start index of observations
while y(1, idx_current_obs_start) == k_first_measurement
    idx_current_obs_start = idx_current_obs_start + 1;
    if(idx_current_obs_start > length(y(1,:)))
        break;
    end
end
idx_prev_obs_end = idx_current_obs_start - 1;
idx_current_obs_end = idx_current_obs_start;

% Animation
if(visualize)
    %x_gt_min = floor(p_k_groundtruth(1,1) )-5;
    %x_gt_max = ceil(p_k_groundtruth(1,1) )+5;
    x_gt_min = -5;
    x_gt_max = 35;
    figure
    hold on
    grid on
    xlim([-5 30]);
    ylim([-2 20]);
    h_p_k = plot(p_k(1,1), 0, 'ro', 'MarkerSize', 5);
    h_p_k_gt = plot(p_k(1,1), 0, 'ko', 'MarkerSize', 10);
    p_k_particles = zeros(n_particles, 1);
    for i = 1:n_particles
        p_k_particles(i) = p_k__i(1, 1, i); 
    end
    h_p_k_particles = plot(p_k_particles, -0.5, 'r.');
    h_birth = plot(0, 0, 'm-');
    h_obs = plot(0, 0, 'b-');
    h_updated = plot(0, 0, 'r-');
    h_obs_lim_max = line( [p_k(1,1) + y_rangeLim  p_k(1,1) + y_rangeLim], [-0.25 0.25] );
    h_obs_lim_min = line( [p_k(1,1) - y_rangeLim  p_k(1,1) - y_rangeLim], [-0.25 0.25] );
    for m = 1:n_features
        line(  [map(1,m) map(1,m)], [-1 0], 'Color', 'k' );
    end
    
    [max_weight i_max_weight] = max(particle_weight(:,1));
end

k_sim_start = k_first_measurement + 1;
k_sim_end = k_max;
for k = k_sim_start:k_sim_end;
    
    if mod(k,50) == 0
        fprintf('Timestep k = %d\n', k);
    end
    if idx_current_obs_start > length(y(1,:))
        idx_current_obs_start = idx_current_obs_start - 1;
        idx_current_obs_end = idx_current_obs_start;
    end
    if y(1, idx_current_obs_start) == k

        while y(1, idx_current_obs_end) < k + 1
            idx_current_obs_end = idx_current_obs_end + 1;
            if(idx_current_obs_end > length(y(1,:)))
                break;
            end
        end
        idx_current_obs_end = idx_current_obs_end - 1;
    
    end
    
    %% Vehicle Pose Prediction (proposal distribution)
    
    t_propogation = tic;
    d_k_km = d_k(:,k);
    D_k_km = reshape(D_vec_k(:,k), 3, 3);

    for i = 1:n_particles
        
        % Pose of particle i at previous timestep k-1
        p_km = p_k__i(:, k-1, i);
        C_km = reshape(C_vec_k__i(:, k-1, i), 3, 3);
        
        % Displacement sampling
        eTran_k_km__i = [noise_motion(1)*randn; noise_motion(2)*randn; noise_motion(3)*randn];
        eRot_vec_k_km__i = [noise_motion(4)*randn; noise_motion(5)*randn; noise_motion(6)*randn];
        eRot_k_km__i = aa2Rmat(eRot_vec_k_km__i);
        
        if i == 1 % Since we are not sampling a large number of particles, just have one sampled at the mean
            eTran_k_km__i = eTran_k_km__i * 0;
            eRot_k_km__i = eye(3);
        end
        
        d_k_km__i = d_k_km + eTran_k_km__i;
        D_k_km__i = eRot_k_km__i*D_k_km;

        % Propogate particle
        [p_k, C_k] = motionModel(p_km, C_km, d_k_km__i, D_k_km__i);
        p_k__i(:, k, i) = p_k;
        C_vec_k__i(:, k, i) = reshape(C_k, 9, 1);
        
        if i == 1 && force_particle_on_true_trajectory == 1 
            p_k__i(:, k, i) = p_k_groundtruth(:, k);
            C_vec_k__i(:, k, i) = c_k_groundtruth(:, k);
        end
        
    end
    
    C_k_dr_vec = reshape(C_k_dr(:,k-1), 3, 3);
    [p_k_dr(:,k), C_k_dr_vec] = motionModel(p_k_dr(:,k-1), C_k_dr_vec, d_k_km, D_k_km);
    C_k_dr(:,k) = reshape(C_k_dr_vec, 9, 1);
    
    time_propogation(k) = toc(t_propogation);
    
    if(visualize)
        delete(h_p_k);
        delete(h_p_k_gt);
        delete(h_obs_lim_max);
        delete(h_obs_lim_min);
        delete(h_p_k_particles);
        h_p_k = plot(p_k__i(1,k,i_max_weight), 0, 'ro', 'MarkerSize', 5);
        h_p_k_gt = plot(p_k_groundtruth(1,k), -0.5, 'ko', 'MarkerSize', 10);
        p_k_particles = zeros(n_particles, 1);
        for i = 1:n_particles
            p_k_particles(i) = p_k__i(1, k, i); 
        end
        h_p_k_particles = plot(p_k_particles, -0.5, 'r.');
        h_obs_lim_max = line( [p_k__i(1,k,i_max_weight) + y_rangeLim  p_k__i(1,k,i_max_weight) + y_rangeLim], [-0.25 0.25] );
        h_obs_lim_min = line( [p_k__i(1,k,i_max_weight) - y_rangeLim  p_k__i(1,k,i_max_weight) - y_rangeLim], [-0.25 0.25] );
        x_gt_min = min(floor(p_k_groundtruth(1,k) )-5, x_gt_min);
        x_gt_max = max(ceil(p_k_groundtruth(1,k) )+5, x_gt_max);

        pause(0.005);

        if y(1, idx_current_obs_start) ~= k
            if(h_obs) ~= 0
                delete(h_obs);
                h_obs = 0;
            end
        end
    end
    
    if y(1, idx_current_obs_start) == k
    
        %% Map update - Correction for features inside sensing area

        t_update = tic;

        % Keep track of which feature is outside field of view
        %Features_inside_FOV_before_update = cell(i);
        
        % For deciding whether birth Gaussians should be created the next
        % timestep
        % birth_Gaussian_mDist_min = zeros( n_particles, idx_current_obs_end - idx_current_obs_start + 1);
        
        % M_size_before_update = zeros(n_particles, 1);
        % weights_before_update = cell(n_particles, 1);

        for i = 1:n_particles

            % M_size_before_update(i) = M_size(i);
            % weights_before_update{i} = M{i,3}(1:M_size_before_update(i));

            % For every pair of features - measurement, make new Gaussian
            % New Gaussians are appended to the end of M
            % Existing Gaussians will be used to represent missed detections

            p_k = p_k__i(:, k, i);
            C_k = reshape(C_vec_k__i(:, k, i), 3, 3);
            % Features_inside_FOV_before_update{i} = zeros(M_size(i), 1); 

            % new_gaussian_weight_numerator_table = zeros(M_size_before_update(i), idx_current_obs_end - idx_current_obs_start + 1);
            % new_gaussian_weight_table_feature_correspondence = new_gaussian_weight_numerator_table;
            % new_gaussian_mDist_table = inf(M_size_before_update(i), idx_current_obs_end - idx_current_obs_start + 1);
            
            % Find features within measurement range
            r = zeros(M_size_before_update(i), 1); % expected feature range
            r_plus =  y_rangeLim + sensor_limit_upper_buffer;
            r_minus = y_rangeLim - sensor_limit_lower_buffer;
            
            n_features_in_fov = 0;
            idx_features_in_fov = zeros(M_size(i), 1);
            for m = 1:M_size(i) 
                p_m = M{i,1}(:,m);
                p_m_k = measureModel(p_k, C_k, p_m);
                r(m) = norm(p_m_k);
                if(r(m) <= r_plus) 
                    n_features_in_fov = n_features_in_fov + 1;
                    idx_features_in_fov(n_features_in_fov) = m;
                end
            end
            idx_features_in_fov(n_features_in_fov + 1 : end) = [];
            n_measurements = idx_current_obs_end - idx_current_obs_start + 1; 
            measurementLikelihoods = zeros(n_measurements, n_features_in_fov);
            updated_map_pos = cell(n_measurements, n_features_in_fov);
            updated_map_cov = cell(n_measurements, n_features_in_fov);

            for m = 1:n_features_in_fov
                
                % Map feature position covariance
                m_idx = idx_features_in_fov(m);
                p_m = M{i,1}(:,m_idx);
                p_m_k = measureModel(p_k, C_k, p_m);
                P_m_km = M{i,2}(1,m_idx);
                
                for n = 1:n_measurements
                    
                    H = 1;
                    R = noise_obs(1, 1);
                    S = H*P_m_km*H' + R;
                    S_det = det(S); 
                    S_inv = 1/S; 
                    
                    z_idx = idx_current_obs_start + n - 1;
                    p_m_k_act = y(2:4, z_idx); 
                    innov = p_m_k_act - p_m_k;
                    mahalanobis_d = innov' * S_inv * innov;
                    
                    if mahalanobis_d <= 3;
                        measurementLikelihoods(n, m) = (2*pi)^(-0.5) * S_det^(-0.5) * exp(-0.5 * mahalanobis_d); 
                        
                        K = P_m_km*H'/S;
                        u_k = p_m + K*innov;
                        P_m_k = (1 - K*H) * P_m_km;
                        P_m_k = (P_m_k + P_m_k')/2;            
                        %P_m_vec_k = reshape([P_m_k, 0, 0, 0, 0, 0, 0, 0, 0], 9, 1);
                        
                        updated_map_pos{n, m} = u_k;
                        updated_map_cov{n, m} = P_m_k;
                        
                    end
                    
                end
                
            end
            
            % TODO, find all likely pairings of from measurementLikelihoods
            
            % Remove all zero rows
            
            % Create new particle for likely pairs
            % M_size(i) = M_size(i) + 1;
            % M{i,1}(:, M_size(i)) = u_k;
            % M{i,2}(:, M_size(i)) = P_m_vec_k;
               
            % Set particle weight

        end

        map_size_update(k) = M_size(i);
        time_update(k) = toc(t_update);
        %fprintf('Particle %d :: updated map size : %d Gaussians\n', i, M_size(i));
        
        if(visualize)
            x_plot = min(floor(p_k__i(1,k,i_max_weight) - y_rangeLim*1.5), x_gt_min) : 0.05 : max(ceil(p_k__i(1,k,i_max_weight) + y_rangeLim*1.5), x_gt_max);
            y_plot = zeros( 1, length(x_plot) );
            for n = idx_current_obs_start : idx_current_obs_end
                u = p_k__i(1,k,i_max_weight) + y(2, n);
                cov = R;
                y_plot = y_plot + pdf('normal', x_plot, u, sqrt(cov));
            end
            if h_obs ~= 0
                delete(h_obs);
            end 
            h_obs = plot(x_plot, y_plot, 'b-');

            x_plot = x_gt_min : 0.05 : x_gt_max;
            y_plot = zeros( 1, length(x_plot) );
            for m = 1:M_size(i_max_weight)
                u = M{i_max_weight,1}(1,m);
                cov = M{i_max_weight,2}(1,m);
                w = M{i_max_weight,3}(m);
                y_plot = y_plot + w*pdf('normal', x_plot, u, sqrt(cov));
            end
            delete(h_updated);
            h_updated = plot(x_plot, y_plot, 'r-');
            pause(0.005);
        end

        

        %% Particle Resampling

        t_resampling = tic;

        if resample_interval ~= -1

            new_particle_set_correspondence = zeros(n_particles, 1);
            total_weight = 0;
            for i = 1 : n_particles
                total_weight = total_weight + particle_weight(i, k);
            end
            sum_normalized_weight_squared = 0;
            for i = 1 : n_particles
                sum_normalized_weight_squared = sum_normalized_weight_squared + (particle_weight(i, k) / total_weight)^2;
            end
            effective_n_particles = 1 / sum_normalized_weight_squared;

            if(effective_n_particles <= effective_particle_threshold)

                sampling_itvl = total_weight / n_particles;
                sample_offset = rand * sampling_itvl;
                cumulative_weight = 0;
                particles_sampled = 0;
                for i = 1 : n_particles
                    cumulative_weight = cumulative_weight + particle_weight(i, k);
                    while(cumulative_weight >= sample_offset)
                        particles_sampled = particles_sampled + 1;
                        new_particle_set_correspondence(particles_sampled) = i;
                        sample_offset = sample_offset + sampling_itvl;
                        if sample_offset > total_weight
                            break; 
                        end
                    end
                end

                % Check which entries we can overwrite
                particle_set_free_spot_ = ones(n_particles, 1);
                j = 1;
                for i = 1 : n_particles         
                     if i == new_particle_set_correspondence(j)
                         % particle i will not get overwritten
                         particle_set_free_spot_(i) = 0;
                         new_particle_set_correspondence(j) = -1; 
                         j = j + 1;
                         if j > n_particles
                            break;
                         end
                         while i == new_particle_set_correspondence(j)
                             j = j + 1; 
                             if j > n_particles
                                 break;
                             end
                         end 
                         if j > n_particles
                            break;
                         end
                    end     
                end

                %display(particle_weight(:,k));
                %display(new_particle_set_correspondence);
                %display(particle_set_free_spot_);

                % Do the actual copying
                free_spot_idx_ = 1;
                for i = 1 : n_particles 
                    j = new_particle_set_correspondence(i);
                    if j >= 1
                        % put particle j from old set into a free spot in new set
                        while particle_set_free_spot_(free_spot_idx_) == 0
                            free_spot_idx_ = free_spot_idx_ + 1; 
                        end
                        particle_set_free_spot_(free_spot_idx_) = 0;

                        % Copy particle pose
                        p_k__i(:, :, free_spot_idx_) = p_k__i(:, :, j);
                        C_vec_k__i(:, :, free_spot_idx_) = C_vec_k__i(:, :, j);
                        % Copy map
                        M{free_spot_idx_, 1} = M{j, 1};
                        M{free_spot_idx_, 2} = M{j, 2};
                        M{free_spot_idx_, 3} = M{j, 3};
                        M_size(free_spot_idx_) = M_size(j);
                    end
                end

                % Equalize all weights
                for i = 1 : n_particles 
                    particle_weight(i, k) = 1;
                end

            end

        end

        time_resampling(k) = toc(t_resampling);
    end

    %% Calculate Errors

    % Highest weight particle
    highest_weight_particle_idx = 1;
    highest_weight = particle_weight(1,k);
    for i = 2:n_particles
         if particle_weight(i,k) > highest_weight
              highest_weight_particle_idx = i;
              highest_weight = particle_weight(i,k);
         end
    end

    % Trajectory errors

    for i = 1:n_particles
        p_k_weighted(k) = p_k_weighted(k) + particle_weight(i,k) * p_k__i(1, k, i);
    end
    p_k_weighted(k) = p_k_weighted(k) / sum(particle_weight(:,k));
    p_k_max_weight(k) = p_k__i(1, k, highest_weight_particle_idx); 

    % Map Cardinality error

    for i = 1:n_particles
        nFeaturesEstimate(i,k) = sum(M{i,3}(1:M_size(i)));
    end
    nFeaturesEstimateAllParticles(k) = sum (particle_weight(:,k) .* nFeaturesEstimate(:,k)) / sum(particle_weight(:,k));
    nFeaturesEstimateMaxWeightParticle(k) = nFeaturesEstimate(highest_weight_particle_idx,k);
    for n = idx_current_obs_start : idx_current_obs_end
        idx = y(5,n);
        if(idx > 0)
            featuresObserved(idx) = 1;
        end
    end
    nFeaturesObserved(k) = sum(featuresObserved);
    
    % Map Wasserstein-distance error using highest-weight particle
    
    map_errors = zeros(n_particles, 2);
    map_estimate_error_all_particles(k) = 0;
    c2 = cutoff ^ 2;
    for i = 1:n_particles
        features_estimated = 0;
        dist_error2 = 0;
        for j = 1:M_size(i) % estimated map
            e = M{i,1}(:,j) - map(:,1);
            w = round(M{i,3}(j));
            features_estimated = features_estimated + w;
            min_dist2 = e'*e;
            for m = 2:n_features % find closest groundtruth feature
                e = M{i,1}(:,j) - map(:,m); 
                dist2 = e'*e;
                min_dist2 = min(min_dist2, dist2);
            end
            min_dist2 = min(min_dist2, c2) * w;
            dist_error2 = dist_error2 + min_dist2;
        end
        dim_error2 = c2 * abs(nFeaturesObserved(k) - features_estimated);
        map_errors(i,1) = sqrt( (dist_error2 + dim_error2) / nFeaturesObserved(k) );
        map_errors(i,2) = particle_weight(i,k);
        map_estimate_error_all_particles(k) = map_estimate_error_all_particles(k) + map_errors(i,1) * map_errors(i,2);
    end
    map_estimate_error_all_particles(k) = map_estimate_error_all_particles(k) / sum(map_errors(:,2));
    [highest_particle_weight, i] = max(particle_weight(:,k-1));
    map_estimate_error_highest_weight_particle(k) = map_errors(i,1);
    
%     c2 = cutoff ^ 2;
%     dist_error_all_particles2 = 0;
%     dist_error_highest_weight_particle2 = 0;
%     for i = 1:n_particles    
%         dist_error2_all_features_for_particle_i = 0;
%         for j = 1:M_size(i) % estimated map
%             e = M{i,1}(:,j) - map(:,1);
%             min_dist2 = e'*e;
%             for m = 2:n_features % find closest groundtruth feature
%                 e = M{i,1}(:,j) - map(:,m); 
%                 dist2 = e'*e;
%                 min_dist2 = min(min_dist2, dist2);
%             end
%             min_dist2 = min(min_dist2, c2);
%             gaussian_weight = M{i,3}(j);
%             min_dist2 = min_dist2 * gaussian_weight^2;
%             dist_error2_all_features_for_particle_i = dist_error2_all_features_for_particle_i + min_dist2;
%         end
%         dist_error_all_particles2 = dist_error_all_particles2 + particle_weight(i,k)^2 * dist_error2_all_features_for_particle_i;
%         if i == highest_weight_particle_idx
%             dist_error_highest_weight_particle2 = dist_error2_all_features_for_particle_i;
%         end
%     end
%     dist_error_all_particles2 = dist_error_all_particles2 /  sum(particle_weight(:,k))^2;
%     dim_error_all_particles = c2 * abs(nFeaturesObserved(k) - nFeaturesEstimateAllParticles(k));
%     dim_error_highest_weight_particle = c2 * abs(nFeaturesObserved(k) - nFeaturesEstimateMaxWeightParticle(k));
%     map_estimate_error_all_particles(k) = sqrt( (dist_error_all_particles2 + dim_error_all_particles) / nFeaturesObserved(k) );
%     map_estimate_error_highest_weight_particle(k) = sqrt( (dist_error_highest_weight_particle2 + dim_error_highest_weight_particle) / nFeaturesObserved(k) );

    %% Update parameters for next timestep
    idx_prev_obs_start = idx_current_obs_start;
    idx_prev_obs_end = idx_current_obs_end;
    idx_current_obs_start = idx_current_obs_end + 1;
    idx_current_obs_end = idx_current_obs_start;
        
end


