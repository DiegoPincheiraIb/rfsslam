% PHD Filter 
% Keith Leung 2013
%
% Functions to customize:
%   generateTrajectory
%   generateMeasurements
%   motionModel
%   measureModel
%   inverseMeasureModel
%   inverseMeasureUncertainty

% phdFilterSimSetup;

% Vehicle trajectory (6-dof poses) for each particle
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
w_before_update = zeros(n_particles, k_max); % for empty strategy
w_after_update = zeros(n_particles, k_max); % for empty strategy
v_eval_pos = zeros(3, n_particles); % For single feature strategy
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
    x_gt_min = floor(p_k_groundtruth(1,1) )-5;
    x_gt_max = ceil(p_k_groundtruth(1,1) )+5;
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
    
    %% Birth Gaussians (new landmarks)
    % Use measurements from previous timestep to generate birth gaussians
    % The birth set gets added to the predicted map set
    % Since the map is static, predicted map set is the same as the updated
    % map set from the previous timestep
    % Hence we will add the birth gaussians directly to set M
    
    %fprintf('Creating %d birth Gaussians\n', idx_prev_obs_end - idx_prev_obs_start + 1);
    t_birth = tic;
    
    if y(1, idx_prev_obs_start) == k - 1
    
        for i = 1:n_particles

            n_birth = 0;
            for obs_idx = idx_prev_obs_start : idx_prev_obs_end

                % Control birth Gaussians based on likelihood of previous measurements
                
                % We only checked with features that we think were inside
                % the sensor field of view at the previous time step.
                % If a measurement does not associate any any features
                % inside the field of view, we should also check outside.
                
                new_feature_is_likely = 0;
                if isempty(birth_Gaussian_mDist_min) || birth_Gaussian_mDist_min(i, obs_idx - idx_prev_obs_start + 1) > birth_Gaussian_likelihood_threshold
                    new_feature_is_likely = 1;

                    % Check with features that may be outside sensor FOV
                    for m = 1 : M_size(i)
                        
                        p_m_k_act = y(2:4, obs_idx); 
                        p_m_k = M{i,1}(:, m) - p_k__i(:,k-1,1);
                        P_m_k = M{i,2}(1, m);
                        
                        innov = p_m_k_act - p_m_k;
                        H = 1;
                        R = noise_obs(1, 1);
                        S = H*P_m_k*H' + R;

                        g_mahalanobis_dist = innov' / S * innov;
                        
                        if g_mahalanobis_dist < birth_Gaussian_likelihood_threshold
                            new_feature_is_likely = 0;
                            break; % we found a feature that the measurement isw likely to have come from 
                        end
                        
                    end
                    % Without the above in place, many new Gaussians get
                    % created
                end
                
                if new_feature_is_likely == 1

                    M_size(i) = M_size(i) + 1; % Increase the birth set size by 1 for particle i
                    n_birth = n_birth + 1;

                    % Pose of particle i at previous timestep k-1
                    p_km = p_k__i(:, k-1, i);
                    C_km = reshape(C_vec_k__i(:, k-1, i), 3, 3);

                    % A measurement at previous timestep k-1, 1)
                    p_m_km = y(2:4, obs_idx);
                    Sp_m_km = reshape(Y_vec(:, obs_idx), 3, 3);

                    % Birth Gaussian
                    p_m_i = inverseMeasureModel(p_km, C_km, p_m_km); % feature position using inverse measurement model
                    Sp_m_i = inverseMeasureUncertainty(p_km, C_km, Sp_m_km); % feature position uncertainty

                    M{i,1}(:, M_size(i)) = p_m_i;
                    M{i,2}(:, M_size(i)) = reshape(Sp_m_i, 9, 1);
                    M{i,3}(M_size(i)) = 0.25; % constant birth weight - should relate to probability of detection?

                end

            end

            w_before_update(i, k) = sum( M{i,3}(1:M_size(i)) );

        end
        
        if(visualize)
            x_plot = x_gt_min : 0.05 : x_gt_max;
            y_plot = zeros( 1, length(x_plot) );
            for m = 1:M_size(i_max_weight)
                u = M{i_max_weight,1}(1,m);
                cov = M{i_max_weight,2}(1,m);
                w = M{i_max_weight,3}(m);
                y_plot = y_plot + w*pdf('normal', x_plot, u, sqrt(cov));
            end
            delete(h_birth);
            h_birth = plot(x_plot, y_plot, 'm-');
            pause(0.005);
        end
    
    end
    
    time_birth(k) = toc(t_birth);
    map_size_birth(k) = M_size(i);
    %fprintf('%d birth Gaussians created\n', n_birth);
    
    %% Vehicle Pose Prediction (proposal distribution)
    
    %fprintf('Propogating particles\n');
    
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
    end
    
    if(visualize)
        if y(1, idx_current_obs_start) ~= k
            if(h_obs) ~= 0
                delete(h_obs);
                h_obs = 0;
            end
        end
    end
    
    if y(1, idx_current_obs_start) == k
    
        %% Map update - Correction for features inside sensing area

        %fprintf('Updating Map\n');
        t_update = tic;

        % Keep track of which feature is outside field of view
        Features_inside_FOV_before_update = cell(i);
        
        % For deciding whether birth Gaussians should be created the next
        % timestep
        birth_Gaussian_mDist_min = zeros( n_particles, idx_current_obs_end - idx_current_obs_start + 1);
        
        M_size_before_update = zeros(n_particles, 1);
        weights_before_update = cell(n_particles, 1);

        for i = 1:n_particles

            M_size_before_update(i) = M_size(i);
            weights_before_update{i} = M{i,3}(1:M_size_before_update(i));

            % For every pair of features - measurement, make new Gaussian
            % New Gaussians are appended to the end of M
            % Existing Gaussians will be used to represent missed detections

            p_k = p_k__i(:, k, i);
            C_k = reshape(C_vec_k__i(:, k, i), 3, 3);
            Features_inside_FOV_before_update{i} = zeros(M_size(i), 1); 

            new_gaussian_weight_numerator_table = zeros(M_size_before_update(i), idx_current_obs_end - idx_current_obs_start + 1);
            new_gaussian_weight_table_feature_correspondence = new_gaussian_weight_numerator_table;
            new_gaussian_mDist_table = inf(M_size_before_update(i), idx_current_obs_end - idx_current_obs_start + 1);
            
            r = zeros(M_size_before_update(i), 1); % expected feature range
            r_plus =  y_rangeLim + sensor_limit_upper_buffer;
            r_minus = y_rangeLim - sensor_limit_lower_buffer;

            for m = 1:M_size_before_update(i) 

                % Determine probability of detection for feature m    
                % inside sensing range
                p_m = M{i,1}(:,m);
                p_m_k = measureModel(p_k, C_k, p_m);
                r(m) = norm(p_m_k);

                P_detection = 0;
                if(r(m) <= r_plus) 
                    % inside sensing area
                    Features_inside_FOV_before_update{i}(m) = 1;
                    P_detection = P_detection_static;

                end
                P_missed_detection = 1 - P_detection;

                if P_detection > 0 % Do not need to copy and update feature position if P_detection = 0

                    % Map feature position covariance
                    P_m_km = M{i,2}(1,m);

                    % Map feature weight
                    w_km = M{i,3}(m);
                    M{i,3}(m) = w_km * P_missed_detection; % Update weight for missed detection

                    % measurement model Jacobian
                    % expected measurement model is p_m_k = C_k * (p_m - p_k) + noise;  
                    % Note that, p_k is a particle (sample)
                    H = 1; % *Measure model dependent, move inside for-loop if measuremnet dependent*

                    % measurement noise
                    R = noise_obs(1, 1); % *Measure model dependent, move inside for-loop if measurement dependent*

                    % innovation convariance, *move inside for-loop if R is measurement dependent 
                    S = H*P_m_km*H' + R;
                    S_det = det(S); % Will use later for measurement likelihood
                    S_inv = 1/S; % Will use later for measurement likelihood

                    % Kalman gain, *move inside for-loop if R is measurement dependent 
                    K = P_m_km*H'/S;

                    % Updated covariance
                    P_m_k = (1 - K*H) * P_m_km;
                    P_m_k = (P_m_k + P_m_k')/2;            
                    P_m_vec_k = reshape([P_m_k, 0, 0, 0, 0, 0, 0, 0, 0], 9, 1);

                    for n = idx_current_obs_start : idx_current_obs_end
                        % for each measurement of this timestep

                        % We are generating a new Gaussian for each measurement
                        % of feature m, and adding this feature to the end of M

                        % innovation : actual - predicted measurement
                        p_m_k_act = y(2:4, n); 
                        innov = p_m_k_act - p_m_k;

                        % Updated map feature m's mean and covariance
                        u_k = p_m + K*innov;

                        % Record updated estimates
                        M_size(i) = M_size(i) + 1;
                        M{i,1}(:, M_size(i)) = u_k;
                        M{i,2}(:, M_size(i)) = P_m_vec_k;
                        %plot3(u_k(1), u_k(2), u_k(3), 'b.'); 

                        % Determine weight (numerator) - we can't get actual
                        % weight until we get all weights associated with same
                        % measurement, so we have to store numerators in a table
                        g_mahalanobis_dist = innov' * S_inv * innov;
                        g = (2*pi)^(-0.5) * S_det^(-0.5) * exp(-0.5 * g_mahalanobis_dist); % measurement likelihood

                        % Table (m x n) for weighting factor numerator
                        new_gaussian_weight_numerator_table(m, n - idx_current_obs_start + 1) = w_km * P_detection * g;
                        new_gaussian_mDist_table(m, n - idx_current_obs_start + 1) = g_mahalanobis_dist; 

                        % This is the index number of the new feature in set M
                        new_gaussian_weight_table_feature_correspondence(m, n - idx_current_obs_start + 1) = M_size(i);

                    end
                    
                end
            end
            
            % sum each column, i.e, for each measurement, total weights for all features
            weights_total = sum(new_gaussian_weight_numerator_table, 1);
            weight_denom = weights_total + clutter_intensity;

            % Now update the weights
            gaussian_weight_table = new_gaussian_weight_numerator_table * 0;
            for m = 1:M_size_before_update(i)
                for n = 1:idx_current_obs_end - idx_current_obs_start + 1
                    idx = new_gaussian_weight_table_feature_correspondence(m, n);
                    if(idx ~= 0) % idx can be 0 if p_detection = 0, in which case we don't bother creating new Gaussians for this feature
                        M{i,3}(idx) = new_gaussian_weight_numerator_table(m,n) /  weight_denom(n);
                        gaussian_weight_table(m, n) = M{i,3}(idx);
                    end
                end
            end
            
            % Features in the probability of detection ambiguity zone
            for m = 1:M_size_before_update(i)
                if(r(m) <= r_plus && r(m) >= r_minus) 
                    w_minus = weights_before_update{i}(m);
                    w_plus = sum(gaussian_weight_table(m, :)); 
                    w_delta = max( P_detection_static * w_minus - w_plus , 0 ); % Weight taken away by negative information
                    M{i,3}(m) = M{i,3}(m) + w_delta;
                end
            end
            
            w_after_update(i, k) = sum( M{i,3} );   
            
            % For controlling birth Gaussians for the next timestep
            birth_Gaussian_mDist_min(i, :) = min(new_gaussian_mDist_table, [], 1);

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

        
        %% Determine particle weight

        t_weighting = tic;
        for i = 1:n_particles
                 
            % Choose map set for particle
            n_measurements = idx_current_obs_end - idx_current_obs_start + 1; 
            if (particle_weighting_feautre_set_max_size >= 0)
                weighting_map_set_size = min(particle_weighting_feautre_set_max_size, n_measurements);
            else
                weighting_map_set_size = n_measurements;
            end

            if weighting_map_set_size > 0
                
                if particle_weighting_random_map == 0
                
                    weighting_map_set = -ones(weighting_map_set_size, 2); % [map_feature_index weight]
                    smallest_weight_idx = 1;
                    if M_size(i) > M_size_before_update(i)
                        for j = M_size_before_update(i) + 1:M_size(i)
                            w = M{i,3}(j);
                            for m = 1:weighting_map_set_size
                                if w > weighting_map_set(m, 2)
                                    weighting_map_set(smallest_weight_idx, 1) = j;
                                    weighting_map_set(smallest_weight_idx, 2) = w; 
                                    break;
                                end
                            end
                            smallest_weight = weighting_map_set(1, 2);
                            smallest_weight_idx = 1;
                            for m = 2:weighting_map_set_size
                                if weighting_map_set(m, 2) < smallest_weight
                                    smallest_weight = weighting_map_set(m, 2);
                                    smallest_weight_idx = m;
                                end
                            end
                        end
                    end

                    for j = 1:weighting_map_set_size
                        if weighting_map_set(j, 1) < 0.5 % todo move this to setup
                            weighting_map_set(j, 1) = -1;
                            weighting_map_set_size = weighting_map_set_size - 1;
                        end
                    end

                    map_eval_pos = zeros(3, weighting_map_set_size);
                    for j = 1:weighting_map_set_size
                        map_eval_pos(:, j) = M{i,1}(:,weighting_map_set(j, 1));
                    end
                    
                elseif particle_weighting_random_map == 1
                    
                    % Pick random map (to show it doesn't really work)
                    map_eval_pos = zeros(3, weighting_map_set_size);
                    p_k = p_k__i(:, k, i);
                    for j = 1:weighting_map_set_size
                        map_eval_pos(1, j) = map_eval_pos(1, j) + p_k(1) + rand * y_rangeLim * 2 - y_rangeLim; 
                    end
                end

            end
            
            % 1. Feature count difference
            feature_count_factor(i, k) = exp(w_after_update(i, k) - w_before_update(i, k));

            % 2. Map intensity Ratio

            if weighting_map_set_size > 0
                v_before_update(i, k) = 0;
                v_after_update(i, k) = 0;
            else
                v_before_update(i, k) = 1;
                v_after_update(i, k) = 1;
            end

            for j = 1:weighting_map_set_size

                for m = 1:M_size_before_update(i)
                    u = M{i,1}(:,m);
                    S = M{i,2}(1,m);
                    S_inv = 1 / S;
                    w = weights_before_update{i}(m); 
                    d = map_eval_pos(:, j) - u;
                    md = d' * S_inv * d;
                    v_before_update(i, k) = v_before_update(i, k) + w * (2*pi)^(-0.5) * det(S)^(-0.5) * exp(-0.5 * md);
                end

                for m = 1:M_size(i)
                    u = M{i,1}(:,m);
                    S = M{i,2}(1,m);
                    S_inv = 1 / S;
                    w = M{i,3}(m); 
                    d = map_eval_pos(:, j) - u;
                    md = d' * S_inv * d;
                    v_after_update(i, k)  = v_after_update(i, k)  + w * (2*pi)^(-0.5) * det(S)^(-0.5) * exp(-0.5 * md);
                end  

            end

            if v_after_update(i,k) == 0
                v_after_update(i,k) = realmin('single'); 
            end
            if v_before_update(i,k) == 0
                v_before_update(i,k) = realmin('single'); 
            end
            similarity_factor(i,k) = v_before_update(i, k) / v_after_update(i, k);

            % 3. Measurement Likelihood 

            if weighting_map_set_size == 0
                measurement_likelihood_factor(i,k) = clutter_intensity ^ (idx_current_obs_end - idx_current_obs_start + 1) / exp(N_c);
            else
                measurement_likelihood_factor(i,k) = 0;
                p_k = p_k__i(:, k, i);
                C_k = reshape(C_vec_k__i(:, k, i), 3, 3);

                likelihoodTable = zeros(n_measurements, weighting_map_set_size);
                table_feature_idx = 0;
                for j = 1:weighting_map_set_size
                    
                        table_feature_idx = table_feature_idx + 1;

                        p_m = map_eval_pos(:, j);
                        p_m_k = measureModel(p_k, C_k, p_m); % predicted measurement of p_m
                        det_R = det(R(1,1));

                        for n = idx_current_obs_start : idx_current_obs_end 
                            p_m_k_act = y(2:4, n); % actual meaurement
                            d = p_m_k_act - p_m_k;
                            md = d' / R(1,1) * d;
                            g = (2*pi)^(-0.5) * det_R^(-0.5) * exp(-0.5 * md); 
                            likelihoodTable(n - idx_current_obs_start + 1, table_feature_idx) = g;
                        end

                end

                likelihoodThreshold = 0.001; % move this to setup file

                % Look for rows that sum to less than threshold, and
                % consider these measurements as clutter
                measurement_is_clutter = sum(likelihoodTable, 2) >= likelihoodThreshold;
                n_clutter_measurements = 0;
                for j = n_measurements: -1 : 1
                    if sum(likelihoodTable(j,:)) < likelihoodThreshold
                        likelihoodTable(j, :) = [];
                        n_clutter_measurements = n_clutter_measurements + 1;
                        n_measurements = n_measurements - 1;
                    end
                end

                measurement_likelihood_factor(i,k) = P_detection_static^weighting_map_set_size * multiFeatureLikelihood(likelihoodTable, likelihoodThreshold, clutter_intensity, N_c) * clutter_intensity ^ n_clutter_measurements;
                n_missed_detection = 0;
                while(measurement_likelihood_factor(i,k) == 0)
                    if n_measurements <= weighting_map_set_size
                        likelihoodTable = [likelihoodTable ones(n_measurements,1)];
                        n_clutter_measurements = n_clutter_measurements + 1;
                        n_missed_detection = n_missed_detection + 1;
                        measurement_likelihood_factor(i,k) = P_detection_static^(weighting_map_set_size - n_missed_detection) * (1 - P_detection_static)^n_missed_detection * multiFeatureLikelihood(likelihoodTable, likelihoodThreshold, clutter_intensity, N_c)  * clutter_intensity ^ n_clutter_measurements;
                    else
                        likelihoodTable = [likelihoodTable; ones(1, weighting_map_set_size)];
                        n_missed_detection = n_missed_detection + 1;
                        measurement_likelihood_factor(i,k) = P_detection_static^(weighting_map_set_size - n_missed_detection) * (1 - P_detection_static)^n_missed_detection * multiFeatureLikelihood(likelihoodTable, likelihoodThreshold, clutter_intensity, N_c)  * clutter_intensity ^ n_clutter_measurements;
                    end
                end
            end
            

            particle_weight(i, k) = measurement_likelihood_factor(i,k) * similarity_factor(i,k) * feature_count_factor(i, k) * particle_weight(i, k-1);


        end
        
        
        % Scale the weights so that they are close to 1
        particle_weight(:, k) = particle_weight(:, k) / mean(particle_weight(:, k));
        [max_weight i_max_weight] = max(particle_weight(:,k));
        
        time_weighting(k) = toc(t_weighting);


        %% Merging and Pruning of Gaussians
        % Processed Gaussians will be marked with sorted_weights[m] = 0
        % Todo - features outside FOV should not even be checked.
        % Todo - determine cutoff point for weights that are too small
        %        i.e., lower portion of weight_order

        %fprintf('Feature merging and pruning\n')
        t_merging = tic;

        for i = 1:n_particles;

            d_threshold = merging_mahalanoblis_distance_threshold ; % Mahalanobis distance threshold for Gaussian merging
            w_threshold = feature_pruning_weight_threshold;

            % Pre-allocate storage for merged Gaussians
            % This will become the new map set at the end of current timestep
            % Assume worst-case where no Gaussians are merged
            x_merged = zeros(3, map_size_limit);
            S_merged = zeros(9, map_size_limit);
            w_merged = zeros(1, map_size_limit);
            merged_count = 0;
            prune_count_ = 0;

            % PHD Filter 2.0
            % If feature is outside sensor FOV, do not check merging or pruning
            % and automatically move this set to the new merged set
            for m = 1:length(Features_inside_FOV_before_update{i})
                if Features_inside_FOV_before_update{i}(m) == 0
                    merged_count = merged_count + 1;
                    w_merged(:, merged_count) = M{i,3}(m);
                    x_merged(:, merged_count) = M{i,1}(:,m);
                    S_merged(:, merged_count) = M{i,2}(:,m);

                    % Set weight to zero for sorting
                    % Feature will go to end up at the end of the sorted list
                    % for merging and pruning, and will be discarded
                    M{i,3}(m) = 0; 
                end
            end

            % Sort features according to weight
            % Features outside FOV already have weights of 0
            % Features below a certain threshold will get pruned
            
            [sorted_weights, weight_order] = sort(M{i,3}(1:M_size(i)), 2, 'descend');
            for m = 1:length(weight_order)
                if sorted_weights(m) < 0.01
                    sorted_weights(m:end) = [];
                    weight_order(m:end) = [];
                    break;
                end
            end

            % Pre-allocate temporary storage for merging calculations
            % Assume worst-case where all Gaussians are merged into one
            x_ = zeros(3, length(weight_order));
            S_ = zeros(9, length(weight_order));
            w_ = zeros(1, length(weight_order));

            for m = 1:length(weight_order)
                if(sorted_weights(m) ~= 0)  % Do not process if Gaussian has already been used

                    sorted_weights(m) = 0;  % flag that this Gaussian has been used

                    idx_m = weight_order(m);  % Find the index before sorting was performed
                    x_m = M{i,1}(:,idx_m);
                    S_vec_m = M{i,2}(:,idx_m);
                    S_m = S_vec_m(1,1);
                    w_m = M{i,3}(idx_m);

                    gaussians_to_merge_ = 1;
                    x_(:, gaussians_to_merge_) = x_m;
                    S_(:, gaussians_to_merge_) = S_vec_m;
                    w_(:, gaussians_to_merge_) = w_m;
                    
                    for n = 1:length(weight_order)
                        if(sorted_weights(n) ~= 0)  % Do not process if Gaussian has already been used

                            idx_n = weight_order(n);  % Find the index before sorting was performed
                            x_n = M{i,1}(:,idx_n);
                            S_vec_n = M{i,2}(:,idx_n);
                            S_n = S_vec_n(1,1);
                            w_n = M{i,3}(idx_n);

                            % Check if n should be merged with m
                            % If so, store information for processing
                            e = x_m - x_n;
                            dm = e'/S_m*e;
                            dn = e'/S_n*e;

                            if(min(dm, dn) < d_threshold && w_n < merging_weight_threshold)

                                sorted_weights(n) = 0; % flag that this Gaussian has been used

                                gaussians_to_merge_ = gaussians_to_merge_ + 1;
                                x_(:, gaussians_to_merge_) = x_n;
                                S_(:, gaussians_to_merge_) = S_vec_n;
                                w_(:, gaussians_to_merge_) = w_n;

                            end
                        end
                    end

                    % Now do the actual merging and pruning

                    w_merged_ = sum(w_(1:gaussians_to_merge_));

                    if( w_merged_ >= w_threshold)

                        x_merged_ = zeros(3,1);
                        for n = 1:gaussians_to_merge_
                            x_merged_ = x_merged_ + w_(:,n) * x_(:,n);
                        end
                        x_merged_ = x_merged_ / w_merged_;

                        S_merged_ = 0;
                        for n = 1:gaussians_to_merge_
                            w_n = w_(:, n);
                            x_n = x_(:, n);
                            S_n = S_(1, n);
                            d_n = x_merged_ - x_n; d_n = d_n(1);
                            S_merged_ = S_merged_ + w_n*(S_n + feature_merging_covarance_inflation_factor*(d_n*d_n'));
                        end
                        S_merged_ = S_merged_ / w_merged_;
                        S_vec_merged_ = reshape([S_merged_, 0, 0, 0, 0, 0, 0, 0, 0], 9, 1); 
                        

                        merged_count = merged_count + 1;
                        w_merged(:, merged_count) = w_merged_;
                        x_merged(:, merged_count) = x_merged_;
                        S_merged(:, merged_count) = S_vec_merged_;

                        

                    else
                        prune_count_ = prune_count_ + 1;
                    end
                end
            end

            % Copy to overwrite existing data in M

            M_size(i) = merged_count;
            M{i,1} = x_merged;
            M{i,2} = S_merged;
            M{i,3} = w_merged;

        end

        time_merging(k) = toc(t_merging);
        map_size_merging(k) = M_size(i);
        %fprintf('Particle %d :: pruned : %d Gaussians\n', i, prune_count_);
        %fprintf('Particle %d :: new map size : %d Gaussians\n', i, merged_count);



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


