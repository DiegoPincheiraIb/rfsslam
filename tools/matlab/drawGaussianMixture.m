function [gm, h] = drawGaussianMixture(mean, cov, weight, drawLimits, res)

% drawGaussianMixture(mean, cov, weight, drawLimits, res)
%
% Function for drawing Bi-variate Gaussian Mixture
%
% mean - n x 2 array of [x, y] coordinates
% cov  - n x 3 array of [var_xx, var_yy, cov_xy] for covariance 
% weight - n x 1 array of weighting
% drawLime - [x_min, x_max, y_min, x_min]
% res - resolution 

x_min = drawLimits(1);
x_max = drawLimits(2);
y_min = drawLimits(3);
y_max = drawLimits(4);  

n = length(mean(:,1));

gm = zeros(round(x_max - x_min)/res + 1, round(y_max - y_min)/res + 1);

for x = x_min : res : x_max
    idx_x = round((x - x_min)/res + 1); 
    for y = y_min : res : y_max
        idx_y = round((y - y_min)/res + 1);
        for g = 1:n
            S = [cov(g,1), cov(g,3); cov(g,3), cov(g,2)];
            d = [x; y] - (mean(g, :))';
            mahalanobis_d = d' / S * d;
            pdf = 1/((2 * pi) * (det(S))^0.5) * exp(-0.5 * mahalanobis_d); 
            gm(idx_x, idx_y) = gm(idx_x, idx_y) + weight(g) * pdf;
        end
    end
end

h = surf(x_min : res : x_max, y_min : res : y_max, gm);
xLim([x_min, x_max]);
yLim([y_min, y_max]);

 