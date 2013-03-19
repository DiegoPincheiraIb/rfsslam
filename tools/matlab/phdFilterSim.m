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

clear all;
close all;
rand('seed', 1);

% Simulation settings
k_max = 1000;
n_features = 100;
y_rangeLim = 5;
%noise_motion = [0.05, 0.05, 0.05, 0.01, 0.01, 0.01];
noise_motion = 0.000001*[0.005, 0.005, 0.005, 0.001, 0.001, 0.001];
noise_obs = [0.05, 0.01, 0.01];
[p_k, c_k, d_k, D_vec_k] = generateTrajectory(k_max, 0.1, [1, 0.1], noise_motion);
[y, Y_vec, map] = generateMeasurements(p_k, c_k, n_features, y_rangeLim, noise_obs);

figure;
plot3(p_k(1,:), p_k(2,:), p_k(3,:), 'b-');
hold on
plot3(map(1,:), map(2,:), map(3,:), 'k.');
grid on
axis equal

% Filter settings
n_particles = 1;
map_size_limit = 10000; % number of Gaussians in each particle's map

% Vehicle trajectory (6-dof poses) for each particle
p_k__i = zeros(3, k_max, n_particles); 
C_vec_k__i = zeros(9, k_max, n_particles); 
for i = 1:n_particles
    C_vec_k__i(:, 1, i) = reshape(eye(3), 9, 1); 
end

% Per particle map
M = cell(n_particles, 3); % M{:,1} mean, M{:,2} cov reshaped as 9x1 vector, M{:,3} weight 
M_size = zeros(n_particles, 1); % for keeping track of map size
for i = 1:n_particles
    M{i,1} = zeros(3, map_size_limit); % pre-allocate memory for map limit
    M{i,2} = zeros(9, map_size_limit); 
    M{i,3} = zeros(1, map_size_limit);
end

% Performance evaluation
time_birth = zeros(0, k_max);
time_propogation = zeros(0, k_max);
time_update = zeros(0, k_max);
time_sorting = zeros(0, k_max);
time_merging = zeros(0, k_max);
map_size_birth = zeros(0, k_max);
map_size_update = zeros(0, k_max);
map_size_merging = zeros(0, k_max);

%% Some book-keeping
birth_Gaussian_mDist_min = []; % This is used to reduce number of birth Gaussians

k_first_measurement = y(1,1);
idx_prev_obs_start= 1;
idx_current_obs_start= 1;
% find the start index of observations
while y(1, idx_current_obs_start) < k_first_measurement
    idx_current_obs_start = idx_current_obs_start + 1;
    if(idx_current_obs_start > length(y(1,:)))
        break;
    end
end
idx_prev_obs_end = idx_current_obs_start - 1;
idx_current_obs_end = idx_current_obs_start;

k_sim_start = k_first_measurement + 1;
k_sim_end = k_first_measurement + 500;
for k = k_sim_start:k_max % k_sim_end % k_max
    
    fprintf('\nTimestep k = %d\n', k);

    while y(1, idx_current_obs_end) < k + 1
        idx_current_obs_end = idx_current_obs_end + 1;
        if(idx_current_obs_end > length(y(1,:)))
            break;
        end
    end
    idx_current_obs_end = idx_current_obs_end - 1;
    
    %% Birth Gaussians (new landmarks)
    % Use measurements from previous timestep to generate birth gaussians
    % The birth set gets added to the predicted map set
    % Since the map is static, predicted map set is the same as the updated
    % map set from the previous timestep
    % Hence we will add the birth gaussians directly to set M
    
    fprintf('Creating %d birth Gaussians\n', idx_prev_obs_end - idx_prev_obs_start + 1);
    t_birth = tic;
    
    for i = 1:n_particles
        
        n_birth = 0;
        for obs_idx = idx_prev_obs_start : idx_prev_obs_end
            
            % Control birth Gaussians based on likelihood of previous
            % measurements
            if isempty(birth_Gaussian_mDist_min) || birth_Gaussian_mDist_min(obs_idx - idx_prev_obs_start + 1) > 2
            
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
                M{i,3}(M_size(i)) = 1; % constant birth weight - should relate to probability of detection?


                % Plot birth position
                %plot3(p_m_i(1), p_m_i(2), p_m_i(3), 'r.'); 
                % Plot real position
                %c_m = y(5, obs_idx);
                %p_m_i = map(:, c_m);
                %plot3(p_m_i(1), p_m_i(2), p_m_i(3), 'k.'); 
            
            end
            
        end
    end
    
    time_birth(k) = toc(t_birth);
    map_size_birth(k) = M_size(i);
    fprintf('%d birth Gaussians created\n', n_birth);
    
    %% Vehicle Pose Prediction (proposal distribution)
    
    fprintf('Propogating particles\n');
    
    t_propogation = tic;
    d_k_km = d_k(:,k);
    D_k_km = reshape(D_vec_k(:,k), 3, 3);

    for i = 1:n_particles
        
        % Pose of particle i at previous timestep k-1
        p_km = p_k__i(:, k-1, i);
        C_km = reshape(C_vec_k__i(:, k-1, i), 3, 3);
        
        % Displacement sampling
        eTran_k_km__i = [noise_motion(1)*randn; noise_motion(2)*randn; noise_motion(3)*randn];
        d_k_km__i = d_k_km + eTran_k_km__i;
        eRot_k_km__i = [noise_motion(4)*randn; noise_motion(5)*randn; noise_motion(6)*randn];
        D_k_km__i = aa2Rmat(eRot_k_km__i)*D_k_km;

        % Propogate particle
        [p_k, C_k] = motionModel(p_km, C_km, d_k_km__i, D_k_km__i);
        p_k__i(:, k, i) = p_k;
        C_vec_k__i(:, k, i) = reshape(C_k, 9, 1);
    end
    time_propogation(k) = toc(t_propogation);
    
    %% Map update - Correction for features inside sensing area
    
    fprintf('Updating Map\n');
    t_update = tic;
    
    % Keep track of which feature is outside field of view
    Features_inside_FOV_before_update = cell(i);
    
    for i = 1:n_particles

        fprintf('Particle %d :: starting map size : %d Gaussians\n', i, M_size(i));
        M_size_before_update = M_size(i);
        
        % For every pair of features - measurement, make new Gaussian
        % New Gaussians are appended to the end of M
        % Existing Gaussians will be used to represent missed detections
        
        p_k = p_k__i(:, k, i);
        C_k = reshape(C_vec_k__i(:, k, i), 3, 3);
        Features_inside_FOV_before_update{i} = zeros(M_size(i), 1); 
        
        new_gaussian_weight_numerator_table = zeros(M_size_before_update, idx_current_obs_end - idx_current_obs_start + 1);
        new_gaussian_weight_table_feature_correspondence = new_gaussian_weight_numerator_table;
        new_gaussian_mDist_table = inf(M_size_before_update, idx_current_obs_end - idx_current_obs_start + 1);
        
        for m = 1:M_size_before_update 
            
            % Determine probability of detection for feature m    
            % inside sensing range
            p_m = M{i,1}(:,m);
            p_m_k = measureModel(p_k, C_k, p_m);
            r = norm(p_m_k);

            P_detection = 0;
            if(r <= y_rangeLim)      
                % inside sensing area
                Features_inside_FOV_before_update{i}(m) = 1;
                P_detection = 0.95;
            end
            P_missed_detection = 1 - P_detection;
            
            if P_detection > 0 % Do not need to copy and update feature position if P_detection = 0

                % Map feature position covariance
                P_m_km = reshape(M{i,2}(:,m), 3, 3);

                % Map feature weight
                w_km = M{i,3}(m);
                M{i,3}(m) = w_km * P_missed_detection; % Update weight for missed detection

                % measurement model Jacobian
                % expected measurement model is p_m_k = C_k * (p_m - p_k) + noise;  
                % Note that, p_k is a particle (sample)
                H = C_k; % *Measure model dependent, move inside for-loop if measuremnet dependent*

                % measurement noise
                R = diag(noise_obs); % *Measure model dependent, move inside for-loop if measurement dependent*

                % innovation convariance, *move inside for-loop if R is measurement dependent 
                S = H*P_m_km*H' + R;
                S_det = det(S); % Will use later for measurement likelihood
                S_inv = eye(3)/S; % Will use later for measurement likelihood
                

                % Kalman gain, *move inside for-loop if R is measurement dependent 
                K = P_m_km*H'/S;
                
                % Updated covariance
                P_m_k = (eye(3) - K*H) * P_m_km;
                P_m_k = (P_m_k + P_m_k')/2;
                P_m_vec_k = reshape(P_m_k, 9, 1);
                
                % For debug, check P_k is invertible
                P_m_k_inv = inv(P_m_k);

                for n = idx_current_obs_start : idx_current_obs_end % for each measurement of this timestep
                    
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
                    g = (2*pi)^(-3/2) * S_det^(-0.5) * exp(-0.5 * g_mahalanobis_dist); % measurement likelihood
                    
                    % Table (m x n) for weighting factor numerator
                    new_gaussian_weight_numerator_table(m, n - idx_current_obs_start + 1) = w_km * P_detection * g;
                    new_gaussian_mDist_table(m, n - idx_current_obs_start + 1) = g_mahalanobis_dist; 
                    
                    % This is the index number of the new feature in set M
                    new_gaussian_weight_table_feature_correspondence(m, n - idx_current_obs_start + 1) = M_size(i);

                end
                
            end
        end
        
        P_false = 0.02; % update this later
        clutterPHD = P_false; % update this later
        
        % sum each column, i.e, for each measurement, total weights for all features
        weights_total = sum(new_gaussian_weight_numerator_table, 1); 
        birth_Gaussian_mDist_min = min(new_gaussian_mDist_table);
        weight_denom = weights_total + clutterPHD;
        
        % Now update the weights
        for m = 1:M_size_before_update
            for n = 1:idx_current_obs_end - idx_current_obs_start + 1
                idx = new_gaussian_weight_table_feature_correspondence(m, n);
                if(idx ~= 0) % idx can be 0 if p_detection = 0, in which case we don't bother creating new Gaussians for this feature
                    M{i,3}(idx) = new_gaussian_weight_numerator_table(m,n) /  weight_denom(n);
                end
            end
        end
             
    end
    map_size_update(k) = M_size(i);
    time_update(k) = toc(t_update);
    fprintf('Particle %d :: updated map size : %d Gaussians\n', i, M_size(i));
    
    %% Merging and Pruning of Gaussians
    % Processed Gaussians will be marked with sorted_weights[m] = 0
    % Todo - features outside FOV should not even be checked.
    % Todo - determine cutoff point for weights that are too small
    %        i.e., lower portion of weight_order
    
    fprintf('Feature merging and pruning\n')
    t_merging = tic;
    
    for i = 1:n_particles;
        
        d_threshold = 1; % Mahalanobis distance threshold for Gaussian merging
        w_threshold = 0.1; % Weight threshold for pruning
        
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
        % Features below a certain threshold will get pruned, along with
        % features outside the sensor FOV
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
                S_m = reshape(S_vec_m, 3, 3);
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
                        S_n = reshape(S_vec_n, 3, 3);
                        w_n = M{i,3}(idx_n);
                        
                        % Check if n should be merged with m
                        % If so, store information for processing
                        e = x_m - x_n;
                        dm = e'/S_m*e;
                        dn = e'/S_n*e;
                        
                        if(min(dm, dn) < d_threshold)
                            
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
                
                if( w_merged_ > w_threshold)
                    
                    x_merged_ = zeros(3,1);
                    for n = 1:gaussians_to_merge_
                        x_merged_ = x_merged_ + w_(:,n) * x_(:,n);
                    end
                    x_merged_ = x_merged_ / w_merged_;

                    S_merged_ = zeros(3,3);
                    for n = 1:gaussians_to_merge_
                        w_n = w_(:, n);
                        x_n = x_(:, n);
                        S_n = reshape(S_(:, n), 3, 3);
                        d_n = x_merged_ - x_n;
                        S_merged_ = S_merged_ + w_n*(S_n + d_n*d_n');
                    end
                    S_merged_ = S_merged_ / w_merged_;
                    S_vec_merged_ = reshape(S_merged_, 9, 1); 
                    S_merged_inv_ = inv(S_merged_);

                    merged_count = merged_count + 1;
                    w_merged(:, merged_count) = w_merged_;
                    x_merged(:, merged_count) = x_merged_;
                    S_merged(:, merged_count) = S_vec_merged_;
                    
                    %plot3(x_merged(1, merged_count), x_merged(2, merged_count), x_merged(3, merged_count), 'g.'); 
                    
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
    fprintf('Particle %d :: pruned : %d Gaussians\n', i, prune_count_);
    fprintf('Particle %d :: new map size : %d Gaussians\n', i, merged_count);
    
    %% Determine particle weight and resmaple
    for i = 1:n_particles
        
    end
    
    %% Particle Resampling
    
    
    %% Update parameters for next timestep
    idx_prev_obs_start = idx_current_obs_start;
    idx_prev_obs_end = idx_current_obs_end;
    idx_current_obs_start = idx_current_obs_end + 1;
end

%figure;
for i = 1:n_particles
    plot3(p_k__i(1, k_sim_start:k_sim_end, i), p_k__i(2, k_sim_start:k_sim_end, i), p_k__i(3, k_sim_start:k_sim_end, i), 'g-' );
    hold on
end
for m = 1:M_size(1)
    x = M{1,1}(:,m);
    if M{1,3}(m) > 0.5
        plot3(x(1), x(2), x(3), 'ro'); 
    end
end
grid on
axis equal

figure;
hold on
plot(1:k_max, time_birth, 'r');
plot(1:k_max, time_propogation, 'g');
plot(1:k_max, time_update, 'b');
plot(1:k_max, time_merging, 'k');
grid on

figure;
hold on
%plot(1:k_max, map_size_birth, 'r');
plot(1:k_max, map_size_update, 'b');
plot(1:k_max, map_size_merging, 'k');
grid on


    

