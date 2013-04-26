function [y, S_vec] = generate1dMissedDetections(y, S_vec, P_d)

% Function generate1dMissedDetections
% Keith Leung 2013

% Input: y [5xn_measurements] - measurements [timestep, x, y, z, correspondence]
%        S [9xn_measurements] - measurement uncertainty covariance reshaped as
%                        9x1 vector
%        P_d - probability of detection

% Outputs:
%        y [5xn_measurements] - measurements [timestep, x, y, z, correspondence]
%        S [9xn_measurements] - measurement uncertainty covariance reshaped as
%                        9x1 vector

n_measurements = length(y(1,:));
P_md = 1 - P_d;
for m = n_measurements : 1
    if rand() <= P_md
        y(:,m) = [];
        S_vec(:,m) = [];
    end
end
