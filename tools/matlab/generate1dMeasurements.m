function [y, S_vec] = generate1dMeasurements(p_k, c_k, m, limit, noise, Pd, Nc)

% function generateMeasurements
% Keith Leung 2013
%
% Generate feauture measurements along x-axis given a vehicle trajectory
%
% Inputs:
% p_k [3xk] - vehicle positions
% c_k [9xk] - vehicle orientation rotation matrices reshaped to 9x1
% m [3xn_features] - feature positions
% limit - measurement limit parameter, if 0, there is no range limit
% noise - 3x1 uncorrelated measurement additive gaussian noise standard
%         deviations
% Pd - probability of detection
% Nc - Number of clutter expected each timestep (integral of clutter
%      intensity
%
% Outputs:
% y [5xn_measurements] - measurements [timestep, x, y, z, correspondence]
% S [9xn_measurements] - measurement uncertainty covariance reshaped as
%                        9x1 vector

k_max = length(p_k);

% Generate clutter table
nc = 0;
poisson_cmf = zeros(ceil(Nc) + 100, 1);
poisson_cmf(1) = pdf('poiss', 0, Nc);
while poisson_cmf < 0.9995
    nc = nc + 1;
    poisson_pmf = pdf('poiss', nc, Nc);
    poisson_cmf(nc + 1) = poisson_cmf(nc) + poisson_pmf;
end
poisson_cmf(nc : end) = [];


% Generate measurements
n_features = length(m(1,:));
y = zeros(5, (n_features + Nc) *length(p_k)); % preallocate lots of memory
S_vec = zeros(9, (n_features + Nc) * length(p_k));
y_n = 0;
for k = 1:k_max
    for n = 1:n_features
        C_k = reshape(c_k(:,k), 3, 3);
        p_m_k = measureModel(p_k(:,k), C_k, m(:,n));
        measureCov = diag(noise);
        
        % apply measurement limits if any
        if( norm(p_m_k) <= limit || limit == 0)
            
            if( rand <= Pd )
            
                y_n = y_n + 1;
                e = [noise(1)*randn; 0; 0];
                y(:, y_n) = [k; p_m_k + e; n];
                S_vec(:, y_n) = reshape(measureCov, 9, 1);
            
            end
        end
 
    end
    
    % Generate clutter
    randomClutterNum = rand();
    n_clutter_to_generate = 0;
    while randomClutterNum > poisson_cmf(n_clutter_to_generate + 1);
        n_clutter_to_generate = n_clutter_to_generate + 1;
        if(n_clutter_to_generate >= length(poisson_cmf))
            break; 
        end
    end
    for clutter_num = 1:n_clutter_to_generate
        y_n = y_n + 1;
        r = limit;
        while r >= limit
            y_false = [rand()*2*limit - limit; 0; 0];
            r = norm(y_false);
        end
        y(:, y_n) = [k; y_false; -1];
        S_vec(:, y_n) = reshape(measureCov, 9, 1);
    end
    
end

y(:,y_n+1:end) = [];
S_vec(:,y_n+1:end) = [];
