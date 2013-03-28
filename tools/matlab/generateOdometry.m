function [d_k_km_noisy, D_vec_k_km_noisy] = generateOdometry(d_k_km, D_vec_k_km, noise)

% Function generateOdometry
% Keith Leung 2013
%
% Takes displacement inputs and adds noise to them
%
% Inputs:
% d_k_km - noise free vehicle odometry for linear motion
% D_vec_k_km - noise free vehicle odometry for rotation (rotation matrix reshaped to 9x1)
% noise - 6x1 vector of independent motion noise standard deviation
%
% Outputs
% d_k_km_noisy - vehicle odometry for linear motion
% D_vec_k_km_noisy - vehicle odometry for rotation (rotation matrix reshaped to 9x1)

d_k_km_noisy = d_k_km * 0;
D_vec_k_km_noisy = D_vec_k_km * 0;
D_vec_k_km_noisy(:, 1) = reshape(eye(3), 9, 1);

for k = 2 : length(d_k_km(1, :))
   
    noise_translate = [noise(1)*randn; noise(2)*randn; noise(3)*randn];
    d_k_km_noisy(:,k) = d_k_km(:,k) + noise_translate;
    noise_rotation = [noise(4)*randn; noise(5)*randn; noise(6)*randn];
    D_k_km_noisy = aa2Rmat(noise_rotation) * reshape(D_vec_k_km(:, k), 3, 3);
    D_vec_k_km_noisy(:,k) = reshape(D_k_km_noisy, 9, 1);
    
end