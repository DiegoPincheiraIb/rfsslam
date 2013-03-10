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

% Simulation settings
k_max = 1000;
n_features = 100;
y_rangeLim = 5;
noise_motion = [0.1, 0.1, 0.1, 0.01, 0.01, 0.01];
noise_obs = [0.2, 0.1, 0.1];
[p_k, c_k, d_k, D_vec_k] = generateTrajectory(k_max, 0.1, [1, 0.1], noise_motion);
[y, Y_vec, m] = generateMeasurements(p_k, c_k, n_features, y_rangeLim, noise_obs);

figure;
plot3(p_k(1,:), p_k(2,:), p_k(3,:), 'b-');
hold on
plot3(m(1,:), m(2,:), m(3,:), 'k.');
grid on
axis equal

% Filter settings
n_particles = 5;
map_size_limit = 1000; % number of Gaussians in each particle's map

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

%% Some book-keeping
idx_prev_obs_start= 1;
idx_current_obs_start= 1;
% find the start index of observations
while y(1, idx_current_obs_start) < 2
    idx_current_obs_start = idx_current_obs_start + 1;
    if(idx_current_obs_start > length(y(1,:)))
        break;
    end
end
idx_prev_obs_end = idx_current_obs_start - 1;

for k = 2:10 % k_max

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
    
    for i = 1:n_particles
        
        for obs_idx = idx_prev_obs_start : idx_prev_obs_end
            
            M_size(i) = M_size(i) + 1; % Increase the birth set size by 1 for particle i
            
            % Pose of particle i at previous timestep k-1
            p_km = p_k__i(:, k-1, i);
            C_km = reshape(C_vec_k__i(:, k-1, i), 3, 3);
            
            % A measurement at previous timestep k-1
            p_m_km = y(2:4, obs_idx);
            Sp_m_km = reshape(Y_vec(:, obs_idx), 3, 3);
            
            % Birth Gaussian
            p_m_i = inverseMeasureModel(p_km, C_km, p_m_km); % feature position using inverse measurement model
            Sp_m_i = inverseMeasureUncertainty(p_km, C_km, Sp_m_km); % feature position uncertainty
            
            M{i,1}(:, M_size(i)) = p_m_i;
            M{i,2}(:, M_size(i)) = reshape(Sp_m_i, 9, 1);
            M{i,3}(M_size(i)) = 1; % constant birth weight - should relate to probability of detection?
            
        end
    end
    
    %% Vehicle Pose Prediction (proposal distribution)
    
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
    
    
    %% Map update - Correction for features inside sensing area
    
    P_detection = cell(n_particles, 1); % For storing probability of detection
    
    for i = 1:n_particles

        M_size_before_update = M_size(i);
        
        % For every pair of features - measurement, make new Gaussian
        % New Gaussians are appended to the end of M
        % Existing Gaussians will be used to represent missed detections
        
        p_k = p_k__i(:, k, i);
        C_k = reshape(C_vec_k__i(:, k, i), 3, 3);
        P_d{i} = zeros(M_size(i), 1); % assume everything is outside sensing area at first
        
        new_gaussian_weight_table = zeros(M_size_before_update, idx_current_obs_end - idx_current_obs_start + 1);
        new_gaussian_weight_table_feature_correspondence = new_gaussian_weight_table;
        
        for m = 1:M_size_before_update 
            
            % Determine probability of detection
               
            %inside sensing range
            p_m = M{i,1}(:,m);
            p_m_k = measureModel(p_k, C_k, p_m);
            r = norm(p_m_k);

            if(norm(p_m_k) <= y_rangeLim)      
                % inside sensing area
                P_detection{i}(m) = 0.95;
            end
            P_missed_detection = 1 - P_d{i}(m);
            
            if P_missed_detection < 1 % Do not need to copy and update feature position if P_detection = 0

                % Map feature position covariance
                P_km = reshape(M{i,2}(:,m), 3, 3);

                % Map feature weight
                w_km = M{i,3}(m);
                M{i,3}(m) = w_km * P_missed_detection; % Update weight for missed detection

                % measurement model Jacobian
                % measurement model is p_m_k = C_k * (p_m- p_k) + noise;  
                H = C_k; % *Measure model dependent, move inside for-loop if measuremnet dependent*

                % measurement covariance
                R = 0.1*eye(3,3); % *Measure model dependent, move inside for-loop if measuremnet dependent*

                % innovation convariance, *move inside for-loop if R is measurement dependent 
                S = H*P_km*H' + R;

                % Kalman gain, *move inside for-loop if R is measurement dependent 
                K = P_km*H'/S;

                for n = idx_current_obs_start : idx_current_obs_end % for each measurement of this timestep

                    % innovation : actual - predicted measurement
                    innov = y(2:4, n) - p_m_k;

                    % Updated map feature m's mean and covariance
                    u_k = M{i,1}(:,m) + K*innov;
                    P_k = (eye(3) - K*H) * P_km;
                    P_vec_k = reshape(P_k, 9, 1);

                    % Record updated estimates
                    M_size(i) = M_size(i) + 1;
                    M{i,1}(:, M_size(i)) = u_k;
                    M{i,2}(:, M_size(i)) = P_vec_k;

                    % Determine weight (numerator) - we can't get actual
                    % weight until we get all weights associated with same
                    % feature, so we have to store numerators in a table
                    g = (2*pi)^(-3/2) * det(S)^(-0.5) * exp(-0.5 * innov' / S * innov); % measurement likelihood
                    new_gaussian_weight_table(m, n - idx_current_obs_start + 1) = w_km * Pd * g;
                    new_gaussian_weight_table_feature_correspondence(m, n - idx_current_obs_start + 1) = M_size(i);

                end
                
            end
        end
        
        P_false = 0.1; % update this later
        clutterPHD = P_false; % update this later
        weights_total = sum(new_gaussian_weight_table, 1); % sum each column, i.e, for each measurement, total weights for all features
        % Update the weights for the new Gaussian now that we can calculate
        % the denominator of the weight equation
        for m = 1:M_size_before_update
            for n = 1:idx_current_obs_end - idx_current_obs_start + 1
                Weights(:,n) = Weights(:,n) /  (clutterPHD + Weights_total(n));
            end
        end
        
        % Merging and Pruning and Gaussians
        
        
        % Determine particle weight and resmaple
        
        
    end
    
    %% Update parameters for next timestep
    idx_prev_obs_start = idx_current_obs_start;
    idx_prev_obs_end = idx_current_obs_end;
end

figure;
for i = 1:n_particles
    plot3(p_k__i(1, :, i), p_k__i(2, :, i), p_k__i(3, :, i) );
    hold on
end
grid on
axis equal

return;

   
    
    

