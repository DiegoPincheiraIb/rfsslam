function [y, S_vec] = generate1dMeasurements(p_k, c_k, m, limit, noise)

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
%
% Outputs:
% y [5xn_measurements] - measurements [timestep, x, y, z, correspondence]
% S [9xn_measurements] - measurement uncertainty covariance reshaped as
%                        9x1 vector

if(nargin ~= 5)
    error('generateMeasurements:nargChk', 'generateMeasurements takes 5 inputs only');
end
if(length(p_k(:,1)) ~= 3)
    error('generateMeasurements:sizeChk', 'Input p_k must be of dimension 3xk_max');
end
if(length(c_k(:,1)) ~= 9)
    error('generateMeasurements:sizeChk', 'Input c_k must be of dimension 9xk_max');
end
if(length(m(:,1)) ~= 3)
    error('generateMeasurements:sizeChk', 'Input m must be 3xn_features');
end
if(length(p_k(1,:)) ~= length(c_k(1,:)) )
    error('generateMeasurements:sizeChk', 'Numbers of columns of p_k and c_k do not match');
end
if(length(noise) ~= 3 || ~isvector(noise))
    error('generateMeasurements:sizeChk', 'Input p_k must be a vector of length 3');
end

k_max = length(p_k);

% Generate measurements
n_features = length(m(1,:));
y = zeros(5, n_features*length(p_k)); % preallocate lots of memory
S_vec = zeros(9, n_features*length(p_k));
y_n = 0;
for k = 1:k_max
    for n = 1:n_features
        C_k = reshape(c_k(:,k), 3, 3);
        p_m_k = measureModel(p_k(:,k), C_k, m(:,n));
        measureCov = diag(noise);
        
        % apply measurement limits if any
        if( norm(p_m_k) <= limit || limit == 0)
            y_n = y_n + 1;
            e = [noise(1)*randn; 0; 0];
            y(:, y_n) = [k; p_m_k + e; n];
            S_vec(:, y_n) = reshape(measureCov, 9, 1);
        end
 
    end
end

y(:,y_n+1:end) = [];
S_vec(:,y_n+1:end) = [];
