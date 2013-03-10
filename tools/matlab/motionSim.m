% Motion Simulation
% Keith Leung 2013
%
% For testing motion model related functions

clear all
close all

for particle = 1:100

dt = 0.1;
k_max = 100;
x = zeros(3, k_max);
C_k = eye(3);

v = [1; 0; 0.1];
w = [0; 0; 1];

for k = 2:k_max
    v_noise = [randn*0.1; randn*0.1; randn*0.1];
    w_noise = [randn*0.1; randn*0.1; randn*0.1];
    [d, a] = velocities2Displacement(v + v_noise, w + w_noise, dt);
    D = aa2Rmat(a);
    [x(:,k), C_k] = motionModel(x(:,k-1), C_k, d, D);
end

plot3(x(1,:), x(2,:), x(3,:));
hold on

end

axis equal
grid on
