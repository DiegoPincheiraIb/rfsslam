% MH-FastSLAM Filter 
% Keith Leung 2013


phdFilterSimSetup1d;
n_particles_original = n_particles;

% Vehicle trajectory for each particle
p_k__i = zeros(3, k_max, n_particles_max * 50); 
C_vec_k__i = zeros(9, k_max, n_particles_max * 50); 
for i = 1:n_particles_max
    C_vec_k__i(:, 1, i) = reshape(eye(3), 9, 1); 
end
p_k_dr = zeros(3, k_max);
C_k_dr = zeros(9, k_max);
C_k_dr(:,1) = reshape(eye(3), 9, 1); 
p_k_weighted = zeros(1, k_max);
p_k_max_weight = zeros(1, k_max);

% Per particle map
M = cell(n_particles_max * 50, 3); % M{:,1} mean, M{:,2} cov reshaped as 9x1 vector, M{:,3} log-odds for existence
M_size = zeros(n_particles_max * 50, 1); % for keeping track of map size
for i = 1:n_particles_max * 50
    M{i,1} = zeros(3, map_size_limit); % pre-allocate memory for map limit
    M{i,2} = zeros(9, map_size_limit); 
    M{i,3} = zeros(1, map_size_limit);
end
featuresObserved = zeros(n_features, 1);
nFeaturesObserved = zeros(k_max, 1);
nFeaturesEstimate = zeros(n_particles_max * 50, k_max);
nFeaturesEstimateAllParticles = zeros(k_max, 1);
nFeaturesEstimateMaxWeightParticle = zeros(k_max, 1);
map_estimate_error_all_particles = zeros(k_max, 1);
map_estimate_error_highest_weight_particle = zeros(k_max, 1);

Pd = P_detection_static;
Px = P_feature_existence_prior;
Pfa = P_false_alarm_static;
p_existence_given_measurement = ( Pd*Px + (1-Pd)*Pfa*Px ) / ( Pfa + (1-Pfa)*Pd*Px );
p_existence_given_no_measurement = ( (1-Pd)*Px ) / ( (1-Px) + (1-Pd)*Px);
log_odds_existence_given_measurement = log( p_existence_given_measurement / (1 - p_existence_given_measurement) );
log_odds_existence_given_no_measurement = log( p_existence_given_no_measurement / (1 - p_existence_given_no_measurement) );
log_odds_existence_pruning_threshold = log ( feature_pruning_weight_threshold / (1 - feature_pruning_weight_threshold ) );
log_odds_existence_limit = 25;

% Particle weighting
particle_weight = ones(n_particles_max * 50, k_max);

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
        
        n_particles_before_update = n_particles;

        for i = 1:n_particles_before_update;

            p_k = p_k__i(:, k, i);
            C_k = reshape(C_vec_k__i(:, k, i), 3, 3);

            r = zeros(M_size(i), 1); 
            r_plus =  y_rangeLim + sensor_limit_upper_buffer;
            r_minus = y_rangeLim - sensor_limit_lower_buffer;
            
            n_features_in_fov = 0;
            idx_features_in_fov = zeros(M_size(i), 1);
            features_in_fov_buffer_zone = zeros(M_size(i), 1);
            for m = 1:M_size(i) 
                p_m = M{i,1}(:,m);
                p_m_k = measureModel(p_k, C_k, p_m);
                r(m) = norm(p_m_k);
                if(r(m) <= r_plus) 
                    n_features_in_fov = n_features_in_fov + 1;
                    idx_features_in_fov(n_features_in_fov) = m;
                    if(r(m) >= r_minus)
                        features_in_fov_buffer_zone(m) = 1;
                    end
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
                    
                    if mahalanobis_d <= data_association_gating_mahalanobis_distance;
                        measurementLikelihoods(n, m) = (2*pi)^(-0.5) * S_det^(-0.5) * exp(-0.5 * mahalanobis_d); 
                        
                        K = P_m_km*H'/S;
                        u_k = p_m + K*innov;
                        P_m_k = (1 - K*H) * P_m_km;
                        P_m_k = (P_m_k + P_m_k')/2;            
                        %P_m_vec_k = reshape([P_m_k, 0, 0, 0, 0, 0, 0, 0, 0], 9, 1);
                        
                        updated_map_pos{n, m} = u_k;
                        updated_map_cov{n, m} = eye(3);
                        updated_map_cov{n, m}(1,1) = P_m_k;
                        
                    end
                    
                end
                
            end
        
            
            if n_features_in_fov == 0
                for n = 1:n_measurements
                    z_idx = idx_current_obs_start + n - 1;
                    p_m_k_act = y(2:4, z_idx); 
                    p_m = p_k + p_m_k_act;
                    R = zeros(3); R(1,1) = noise_obs(1);
                    P_m = R;
                    M_size(i) = M_size(i) + 1;
                    M{i,1}(:, M_size(i)) = p_m;
                    M{i,2}(:, M_size(i)) = reshape( P_m, 9, 1 );
                    M{i,3}(M_size(i)) = log_odds_existence_given_measurement; % prob of existence
                end
                sequence = [];
            else
                n_associated_measurement = 0;
                measurementLikelihoods_pruned = zeros( size( measurementLikelihoods ) );
                idx_associated_measurement = zeros(n_measurements, 1);
                for n = 1:n_measurements
                    if sum(measurementLikelihoods(n,:), 2) == 0
                        z_idx = idx_current_obs_start + n - 1;
                        p_m_k_act = y(2:4, z_idx); 
                        p_m = p_k + p_m_k_act;
                        R = zeros(3); R(1,1) = noise_obs(1);
                        P_m = R;
                        M_size(i) = M_size(i) + 1;
                        M{i,1}(:, M_size(i)) = p_m;
                        M{i,2}(:, M_size(i)) = reshape( P_m, 9, 1 );
                        M{i,3}(M_size(i)) = log_odds_existence_given_measurement; 
                    else
                        n_associated_measurement = n_associated_measurement + 1;
                        idx_associated_measurement(n_associated_measurement) = n;
                        measurementLikelihoods_pruned( n_associated_measurement, : ) = measurementLikelihoods(n, :);  
                    end
                end
                measurementLikelihoods_pruned(n_associated_measurement + 1 : end, : ) = [];
                
                % possible associations for each measurement / feature
                n_likely_features = sum(measurementLikelihoods_pruned > 0, 2);
                n_likely_measurements = sum(measurementLikelihoods_pruned > 0, 1);
                
                % initiate data association sequence for checking
                if n_features_in_fov >= n_associated_measurement 
                    sequence = 1:n_features_in_fov;
                else
                    sequence = 1:n_associated_measurement;                
                end
                
                template_particle_i_idx = i;
                template_particle_i_map_pos = M{i,1}; 
                template_particle_i_map_cov = M{i,2}; 
                template_particle_i_map_w = M{i,3}; 
                template_particle_i_map_size = M_size(i); 
                template_particle_i_pos = p_k__i(:, :, i);
                template_particle_i_rot = C_vec_k__i(:, :, i);
                template_particle_i_w = particle_weight(i,k-1);

                % Find top data association sequences
                n_data_association_sequences = 0;
                top_association_sequences = cell(n_particles_max_spawn, 1);
                top_association_likelihoods = zeros(n_particles_max_spawn, 1);
                if n_features_in_fov >= n_associated_measurement
                    while n_data_association_sequences == 0
                        while( ~isempty(sequence) )
                            [sequence, likelihood, sequence_n] = multiFeatureLikelihoodFastSLAM(measurementLikelihoods_pruned, sequence);
                            jointMeasurementLikelihood = prod(likelihood);
                            if jointMeasurementLikelihood > 0
                                n_data_association_sequences = n_data_association_sequences + 1;
                                [min_top_likelihood, min_top_likelihood_idx] = min(top_association_likelihoods);
                                if jointMeasurementLikelihood > min_top_likelihood
                                    top_association_likelihoods( min_top_likelihood_idx ) = jointMeasurementLikelihood;
                                    top_association_sequences{min_top_likelihood_idx} = sequence;
                                end
                            end
                            sequence = sequence_n;       
                        end
                        if  n_data_association_sequences == 0
                            measurementLikelihoods_pruned = [measurementLikelihoods_pruned ones(n_associated_measurement,1)];
                            sequence = 1:length(measurementLikelihoods_pruned(1,:));
                        end
                    end
                else
                    while n_data_association_sequences == 0
                        while( ~isempty(sequence) )
                            [sequence, likelihood, sequence_n] = multiFeatureLikelihoodFastSLAM(measurementLikelihoods_pruned, sequence);
                            jointMeasurementLikelihood = prod(likelihood);
                            if jointMeasurementLikelihood > 0
                                n_data_association_sequences = n_data_association_sequences + 1;
                                [min_top_likelihood, min_top_likelihood_idx] = min(top_association_likelihoods);
                                if jointMeasurementLikelihood > min_top_likelihood
                                    top_association_likelihoods( min_top_likelihood_idx ) = jointMeasurementLikelihood;
                                    top_association_sequences{min_top_likelihood_idx} = sequence;
                                end
                            end
                            sequence = sequence_n;       
                        end
                        if  n_data_association_sequences == 0
                            measurementLikelihoods_pruned = [measurementLikelihoods_pruned; ones(1,n_features_in_fov)];
                            sequence = 1:length(measurementLikelihoods_pruned(:,1));
                        end
                    end
                end
                n_data_association_sequences = min(n_data_association_sequences, n_particles_max_spawn);
                
                % Go through top data association sequences
                if n_features_in_fov >= n_associated_measurement
                    for seq = 1:n_data_association_sequences
                        jointMeasurementLikelihood = top_association_likelihoods(seq);
                        sequence = top_association_sequences{seq};
                        particle_idx = i;
                        if seq > 1
                            % create new particle
                            n_particles = n_particles + 1;
                            particle_idx = n_particles;
                            M{particle_idx,1} = template_particle_i_map_pos; 
                            M{particle_idx,2} = template_particle_i_map_cov; 
                            M{particle_idx,3} = template_particle_i_map_w;
                            M_size(particle_idx) = template_particle_i_map_size;
                            p_k__i(:, :, particle_idx) = template_particle_i_pos; 
                            C_vec_k__i(:, :, particle_idx) = template_particle_i_rot;
                        end
                        particle_weight(particle_idx,k) = jointMeasurementLikelihood * template_particle_i_w;
                        feature_updated = zeros(n_features_in_fov, 1);
                        for n = 1:n_associated_measurement
                            m = sequence(n);
                            if m <= n_features_in_fov
                                % update existing landmark
                                M{particle_idx,1}(:,idx_features_in_fov(m)) = updated_map_pos{idx_associated_measurement(n), m};
                                M{particle_idx,2}(:,idx_features_in_fov(m)) = reshape(updated_map_cov{idx_associated_measurement(n), m}, 9, 1);
                                M{particle_idx,3}(idx_features_in_fov(m)) = M{particle_idx,3}(idx_features_in_fov(m)) + log_odds_existence_given_measurement; 
                                M{particle_idx,3}(idx_features_in_fov(m)) = min( M{particle_idx,3}(idx_features_in_fov(m)), log_odds_existence_limit);
                                feature_updated(m) = 1;
                            else
                                % create new landmark
                                z_idx = idx_current_obs_start + n - 1;
                                p_m_k_act = y(2:4, z_idx); 
                                p_m = p_k + p_m_k_act;
                                P_m = C_k + R;
                                M_size(particle_idx) = M_size(particle_idx) + 1;
                                M{particle_idx,1}(:, M_size(particle_idx)) = p_m;
                                M{particle_idx,2}(:, M_size(particle_idx)) = reshape( P_m, 9, 1 );
                                M{particle_idx,3}(M_size(particle_idx)) = log_odds_existence_given_measurement;
                            end
                        end
                        for m = 1:n_features_in_fov
                            if feature_updated(m) == 0 && features_in_fov_buffer_zone(idx_features_in_fov(m)) ~= 1;
                                 M{particle_idx,3}(idx_features_in_fov(m)) =  M{particle_idx,3}(idx_features_in_fov(m)) + log_odds_existence_given_no_measurement;
                            end
                        end
                        for m = n_features_in_fov : -1 : 1
                            if M{particle_idx,3}(idx_features_in_fov(m)) < log_odds_existence_pruning_threshold;
                                M{particle_idx,1}(:, idx_features_in_fov(m)) = [];
                                M{particle_idx,2}(:, idx_features_in_fov(m)) = [];
                                M{particle_idx,3}(idx_features_in_fov(m)) = []; 
                                M_size(particle_idx) = M_size(particle_idx) - 1;
                            end
                        end
                    end
                else
                    for seq = 1:n_data_association_sequences
                        jointMeasurementLikelihood = top_association_likelihoods(seq);
                        sequence = top_association_sequences{seq};
                        particle_idx = i;
                        if n_data_association_sequences > 1
                            % create new particle
                            n_particles = n_particles + 1;
                            particle_idx = n_particles;
                            M{particle_idx,1} = template_particle_i_map_pos; 
                            M{particle_idx,2} = template_particle_i_map_cov; 
                            M{particle_idx,3} = template_particle_i_map_w;
                            M_size(particle_idx) = template_particle_i_map_size;
                            p_k__i(:, :, particle_idx) = template_particle_i_pos; 
                            C_vec_k__i(:, :, particle_idx) = template_particle_i_rot;
                        end
                        particle_weight(particle_idx,k) = jointMeasurementLikelihood * template_particle_i_w;
                        feature_updated = zeros(n_features_in_fov, 1);
                        used_measurement = zeros(n_associated_measurement, 1);
                        for m = 1:n_features_in_fov
                            n = sequence(m); 
                            if n <= n_associated_measurement
                                used_measurement(n) = 1;
                                % update existing landmark
                                M{particle_idx,1}(:,idx_features_in_fov(m)) = updated_map_pos{idx_associated_measurement(n), m};
                                M{particle_idx,2}(:,idx_features_in_fov(m)) = reshape(updated_map_cov{idx_associated_measurement(n), m}, 9, 1);
                                M{particle_idx,3}(idx_features_in_fov(m)) = M{particle_idx,3}(idx_features_in_fov(m)) + log_odds_existence_given_measurement;
                                M{particle_idx,3}(idx_features_in_fov(m)) = min( M{particle_idx,3}(idx_features_in_fov(m)), log_odds_existence_limit);
                                feature_updated(m) = 1;
                            else
                                % no measurement was associated with m
                                feature_updated(m) = 0;
                            end
                            for n = 1:n_associated_measurement
                                if used_measurement(n) == 0
                                    % create new landmark
                                    z_idx = idx_current_obs_start + n - 1;
                                    p_m_k_act = y(2:4, z_idx); 
                                    p_m = p_k + p_m_k_act;
                                    P_m = C_k + R;
                                    M_size(particle_idx) = M_size(particle_idx) + 1;
                                    M{particle_idx,1}(:, M_size(particle_idx)) = p_m;
                                    M{particle_idx,2}(:, M_size(particle_idx)) = reshape( P_m, 9, 1 );
                                    M{particle_idx,3}(M_size(particle_idx)) = log_odds_existence_given_measurement;
                                end
                            end
                        end
                        for m = 1:n_features_in_fov
                            if feature_updated(m) == 0 && features_in_fov_buffer_zone(idx_features_in_fov(m)) ~= 1
                                 M{particle_idx,3}(idx_features_in_fov(m)) =  M{particle_idx,3}(idx_features_in_fov(m)) + log_odds_existence_given_no_measurement; 
                            end
                        end
                        for m = n_features_in_fov : -1 : 1
                            if M{particle_idx,3}(idx_features_in_fov(m)) < log_odds_existence_pruning_threshold
                                M{particle_idx,1}(:, idx_features_in_fov(m)) = [];
                                M{particle_idx,2}(:, idx_features_in_fov(m)) = [];
                                M{particle_idx,3}(idx_features_in_fov(m)) = []; 
                                M_size(particle_idx) = M_size(particle_idx) - 1;
                            end
                        end 
                    end
                end      
                
                % rescale weights of new particles
                particle_weight(i,k) = particle_weight(i,k) / n_data_association_sequences;
                for j = 0:n_data_association_sequences - 2
                    particle_weight(n_particles - j,k) = particle_weight(n_particles - j,k) / n_data_association_sequences;
                end
                
            end
        end

        map_size_update(k) = M_size(i);
        time_update(k) = toc(t_update);
        %fprintf('Particle %d :: updated map size : %d Gaussians\n', i, M_size(i));
        
        if(visualize)
            x_plot = min(floor(p_k__i(1,k,i_max_weight) - y_rangeLim*1.5), x_gt_min) : 0.05 : max(ceil(p_k__i(1,k,i_max_weight) + y_rangeLim*1.5), x_gt_max);
            y_plot = zeros( 1, length(x_plot) );
            for n = idx_current_obs_start : idx_current_obs_end
                u = p_k__i(1,k,i_max_weight) + y(2, n);
                R = noise_obs(1,1);
                y_plot = y_plot + pdf('normal', x_plot, u, sqrt(R));
            end
            if h_obs ~= 0
                delete(h_obs);
            end 
            h_obs = plot(x_plot, y_plot, 'b-');

            x_plot = x_gt_min : 0.05 : x_gt_max;
            y_plot = zeros( 1, length(x_plot) );
            for m = 1:M_size(i_max_weight)
                u = M{i_max_weight,1}(1,m);
                S = M{i_max_weight,2}(1,m);
                log_w = M{i_max_weight,3}(m);
                w = 1 - 1/(1 + exp(log_w));
                y_plot = y_plot + w*pdf('normal', x_plot, u, sqrt(S));
            end
            delete(h_updated);
            h_updated = plot(x_plot, y_plot, 'r-');
            pause(0.005);
        end

        fprintf('particles: %d\n', n_particles);

        %% Particle Resampling

        t_resampling = tic;

        if resample_interval ~= -1

            new_particle_set_correspondence = zeros(n_particles_original, 1);
            total_weight = 0;
            for i = 1 : n_particles
                total_weight = total_weight + particle_weight(i, k);
            end
            sum_normalized_weight_squared = 0;
            for i = 1 : n_particles
                sum_normalized_weight_squared = sum_normalized_weight_squared + (particle_weight(i, k) / total_weight)^2;
            end
            effective_n_particles = 1 / sum_normalized_weight_squared;
            
            effective_particle_threshold_k = effective_particle_threshold * n_particles / n_particles_original * 2; 
            
            if(effective_n_particles <= effective_particle_threshold_k || n_particles > n_particles_max)
                
                sampling_itvl = total_weight / n_particles_original;
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
                for i = 1 : n_particles_original      
                     if i == new_particle_set_correspondence(j)
                         % particle i will not get overwritten
                         particle_set_free_spot_(i) = 0;
                         new_particle_set_correspondence(j) = -1; 
                         j = j + 1;
                         if j > n_particles_original
                            break;
                         end
                         while i == new_particle_set_correspondence(j)
                             j = j + 1; 
                             if j > n_particles_original
                                 break;
                             end
                         end 
                         if j > n_particles_original
                            break;
                         end
                    end     
                end

                % Do the actual copying
                free_spot_idx_ = 1;
                for i = 1 : n_particles_original
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
                n_particles = n_particles_original;
                for i = 1 : n_particles_original
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
    p_k_weighted(k) = p_k_weighted(k) / sum(particle_weight(1:n_particles,k));
    p_k_max_weight(k) = p_k__i(1, k, highest_weight_particle_idx); 

    % Map Cardinality error

    for i = 1:n_particles
        nFeaturesEstimate(i,k) = 0;
        for m = 1:M_size(i)
            log_w = M{i,3}(m);
            w = 1 - 1/(1 + exp(log_w));
            nFeaturesEstimate(i,k) = nFeaturesEstimate(i,k) + w;
        end
    end
    nFeaturesEstimateAllParticles(k) = sum (particle_weight(1:n_particles,k) .* nFeaturesEstimate(1:n_particles,k)) / sum(particle_weight(1:n_particles,k));
    nFeaturesEstimateMaxWeightParticle(k) = nFeaturesEstimate(highest_weight_particle_idx,k);
    for n = idx_current_obs_start : idx_current_obs_end
        idx = y(5,n);
        if(idx > 0)
            featuresObserved(idx) = 1;
        end
    end
    nFeaturesObserved(k) = sum(featuresObserved);
    
    % Map Wasserstein-distance error 
    
    map_errors = zeros(n_particles, 2);
    map_estimate_error_all_particles(k) = 0;
    c2 = cutoff ^ 2;
    for i = 1:n_particles
        features_estimated = 0;
        dist_error2 = 0;
        for j = 1:M_size(i) % estimated map
            e = M{i,1}(:,j) - map(:,1);
            log_w = M{i,3}(j);
            w = 1 - 1/(1 + exp(log_w));
            w = round(w);
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
    [highest_particle_weight, i] = max(particle_weight(1:n_particles,k-1));
    map_estimate_error_highest_weight_particle(k) = map_errors(i,1);

    %% Update parameters for next timestep
    idx_prev_obs_start = idx_current_obs_start;
    idx_prev_obs_end = idx_current_obs_end;
    idx_current_obs_start = idx_current_obs_end + 1;
    idx_current_obs_end = idx_current_obs_start;
        
end


