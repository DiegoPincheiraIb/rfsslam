function [p_k, c_k, d_k_km, D_vec_k_km] = generateTrajectory(k_max, dt, v_max)

% function generateTrajectory
% Keith Leung 2013
%
% Generate a smooth trajectory in 3d
%
% Inputs:
% k_max - length of trajectory
% dt - size of a timestep
% v_max - [max_linear_velocity max_angular_velocity] per axis
% noise - 6x1 vector of independent motion noise standard deviation
%
% Outputs
% x_k - vehicle positions (3 x k_max)
% p_k - vehicle orientatations (9 x k_max) (rotation matrix reshaped as
%       vectors)
% d_k_km - vehicle odometry for linear motion
% D_vec_k_km - vehicle odometry for rotation (rotation matrix reshaped to
%              9x1)

if(nargin ~= 3)
    error('generateTrajectory:nargChk', 'generateTrajectory takes 3 inputs only');
end
if(~isscalar(k_max))
    error('generateTrajectory:sizeChk', 'Input k_max must be scalar');
end
if(length(v_max) ~= 2 || ~isvector(v_max))
    error('generateTrajectory:sizeChk', 'Input v_max must be a vector of length 2');
end
if(~isscalar(dt))
    error('generateTrajectory:sizeChk', 'Input dt must be scalar');
end

v_max_linear = v_max(1);
v_max_rot = v_max(2);

p_k = zeros(3,k_max);
c_k = zeros(9,k_max);
c_k(:, 1) = [1, 0, 0, 0, 1, 0, 0, 0, 1];
C_k = eye(3);
p_km = zeros(3, 1);
C_km = eye(3);
d_k_km = zeros(3, k_max);
D_vec_k_km = zeros(9, k_max);
D_vec_k_km(:, 1) = [1, 0, 0, 0, 1, 0, 0, 0, 1];


v_k = [(v_max_linear).*rand; 0; 0];
w_k = -v_max_rot + (v_max_rot*2).*rand(3,1);

for k = 2:k_max
    
    if(mod(k, 100) == 0)
       v_k = zeros(3,1);
       while(norm(v_k) < 0.5)
        v_k = [(v_max_linear).*rand; 0; 0];
        w_k = -v_max_rot + (v_max_rot*2).*rand(3,1);
       end
    end
   
    % displacements    
    d_k_km_ = v_k*dt;
    rot_k_km_ = w_k*dt;
    D_k_km_ = aa2Rmat(rot_k_km_);
    
    d_k_km(:, k) = d_k_km_;
    D_vec_k_km(:, k) = reshape(D_k_km_, 9, 1);
    
    % determine pose
    [p_k(:,k) C_k] = motionModel(p_k(:,k-1), C_km, d_k_km_, D_k_km_); 
    c_k(:,k) = reshape(C_k, 9, 1); 
    C_km = C_k;
    
end