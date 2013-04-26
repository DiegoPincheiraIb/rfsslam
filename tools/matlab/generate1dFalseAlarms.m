function[y, S_vec] = generate1dFalseAlarms(y, S_vec, limit, p_fa)

% function generate1dFalseAlarms
% Keith Leung 2013
%
% Input: y [5xn_measurements] - measurements [timestep, x, y, z, correspondence]
%        S [9xn_measurements] - measurement uncertainty covariance reshaped as
%                        9x1 vector
%        limit - measurement limit parameter, if 0, there is no range limit
%        p_fa - probability of false alarm 
%
% Outputs:
%        y [5xn_measurements] - measurements [timestep, x, y, z, correspondence]
%        S [9xn_measurements] - measurement uncertainty covariance reshaped as
%                        9x1 vector

n_measurements = length(y(1,:));
for m = 1:n_measurements
    if rand() <= p_fa
        k = y(1,m);
        r = limit;
        while r >= limit
            y_false = [rand()*2*limit - limit; 0; 0];
            r = norm(y_false);
        end
        y(2:5,m) = [y_false; -1];
    end
end
