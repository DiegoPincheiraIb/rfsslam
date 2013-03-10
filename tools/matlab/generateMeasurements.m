function [y, S_vec, m] = generateMeasurements(p_k, c_k, n_features, limit, noise)

% function generateMeasurements
% Keith Leung 2013
%
% Generate feauture map and measurements given a vehicle trajectory
%
% Inputs:
% p_k [3xk] - vehicle positions
% c_k [9xk] - vehicle orientation rotation matrices reshaped to 9x1
% n_features - number of features to generate along trajectory
% limit - measurement limit parameter
% noise - 3x1 uncorrelated measurement additive gaussian noise standard
%         deviations
%
% Outputs:
% y [4xn_measurements] - measurements [timestep, x, y, z]
% S [9xn_measurements] - measurement uncertainty covariance reshaped as
%                        9x1 vector
% m [3xn_features] - map of observed features

if(nargin ~= 5)
    error('generateMeasurements:nargChk', 'generateMeasurements takes 5 inputs only');
end
if(length(p_k(:,1)) ~= 3)
    error('generateMeasurements:sizeChk', 'Input p_k must be of dimension 3xk_max');
end
if(length(c_k(:,1)) ~= 9)
    error('generateMeasurements:sizeChk', 'Input c_k must be of dimension 9xk_max');
end
if(~isscalar(n_features))
    error('generateMeasurements:sizeChk', 'Input n_features must be scalar');
end
if(length(p_k(1,:)) ~= length(c_k(1,:)) )
    error('generateMeasurements:sizeChk', 'Input p_k must be of dimension 3xk');
end
if(length(noise) ~= 3 || ~isvector(noise))
    error('generateMeasurements:sizeChk', 'Input p_k must be a vector of length 3');
end

k_max = length(p_k);

% Generate landmarks
expectedFeaturesPerStep = n_features / k_max;
n_features = 0;
m = zeros(3, n_features);
for k = 1:k_max
     p_k_i = p_k(:,k);
     C_k_i = reshape(c_k(:,k), 3, 3);
     expectedFeatureCount = expectedFeaturesPerStep * k;
     while(n_features < expectedFeatureCount)
         % generate a measurement
         n_features = n_features + 1;
         p_m_k = -limit + (limit*2).*rand(3,1);
         if(norm(p_m_k) < limit)
             n_features = n_features - 1;
         else
            m(:,n_features) = C_k_i'*p_m_k + p_k_i;
         end
     end
end

% Generate measurements
y = zeros(4, n_features*length(p_k)); % preallocate lots of memory
S_vec = zeros(9, n_features*length(p_k));
y_n = 0;
for k = 1:k_max
    for n = 1:n_features
        C_k = reshape(c_k(:,k), 3, 3);
        p_m_k = measureModel(p_k(:,k), C_k, m(:,n));
        measureCov = [0.2, 0, 0;
                      0, 0.2, 0;
                      0, 0, 0.2];
        
        % apply measurement limits if any
        if( norm(p_m_k) <= limit)
            y_n = y_n + 1;
            e = [noise(1)*randn; noise(2)*randn; noise(3)*randn];
            y(:, y_n) = [k; p_m_k + e];
            S_vec(:, y_n) = reshape(measureCov, 9, 1);
        end
 
    end
end

y(:,y_n+1:end) = [];
S_vec(:,y_n+1:end) = [];
