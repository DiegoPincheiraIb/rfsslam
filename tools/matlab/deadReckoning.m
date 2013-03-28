function [p_k, C_k] = deadReckoning(d_k_km, D_vec_k_km)

% function deadReckoning
% Keith Leung 2013
%
% Function for determining vechile trajectory given inputs (odometry)
%
% Inputs:
% d_k_km - translational displacements (3 x k_max)
% d_vec_k_km - flattened rotation matrices  (9 x k_max)
%
% Outputs:
% p_k - vehicle position (3 x k_max)
% C_k - vechile orientation (flattened rotation matrix 9 x k_max)

p_k = d_k_km * 0;
C_vec_k = D_vec_k_km * 0;
C_vec_k(:, 1) = reshape( eye(3), 9, 1);

k_max = length( d_k_km(1, :) );

for k = 2 : k_max
    C_km = reshape(C_vec_k(:, k-1), 3, 3);
    [p_k(:, k), C_k] = motionModel(p_k(:, k-1), C_km, d_k_km(:, k), reshape(D_vec_k_km(:, k), 3, 3) ); 
    C_vec_k(:, k) = reshape(C_k, 9, 1);
end