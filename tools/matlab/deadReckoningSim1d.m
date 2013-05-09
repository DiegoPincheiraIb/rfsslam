% Dead reckoning 1d monte-carlo simulation
% Keith Leung 2013

figure;
hold on
title('Dead reckoning Monte-Carlo simulation')
plot(1:k_max, p_k_groundtruth(1,:), 'b-', 'LineWidth', 2);
for i = 1:50
    [d_k, D_vec_k] = generateOdometry(d_k_noiseFree, D_vec_k_noiseFree, noise_motion);
    [p_k_dr, C_k_dr] = deadReckoning(d_k, D_vec_k);
    plot(1:k_max, p_k_dr(1, 1:k_max), 'g-' );
end
grid on
