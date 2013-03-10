function [r_k C_k] = motionModel(r_km, C_km, d_k_km, D_k_km)

% Vehicle 3d (Displacement) Motion Model
% Keith Leung 2013
%
% Inputs:
% r_km - 3d position of frame km relative to inertial frame
% C_km - 3d orientation (rotation matrix) of frame km relative to inertial frame 
% d_k_km - 3d displacement from frame km to frame k, expressed in frame km
% C_k_km - 3d rotation (matrix) going from frame km to frame k
%
% Outputs:
% r_k - 3d position of frame k relative to inertial frame
% C_k - 3d orientation (rotation matrix) of frame k relative to inertial frame

if(nargin ~= 4)
    error('motionModel:nargChk', 'aa2Rmat takes 1 input only');
end
if(length(r_km(:,1)) ~= 3 || length(r_km(1,:)) ~= 1)
    error('motionModel:sizeChk', 'Input vector r_km must be 3x1');
end
if(length(C_km(:,1)) ~= 3 || length(C_km(1,:)) ~= 3)
    error('motionModel:sizeChk', 'Input matrix C_km must be 3x3');
end
if(length(d_k_km(:,1)) ~= 3 || length(d_k_km(1,:)) ~= 1)
    error('motionModel:sizeChk', 'Input vector d_k_km must be 3x1');
end
if(length(D_k_km(:,1)) ~= 3 || length(D_k_km(1,:)) ~= 3)
    error('motionModel:sizeChk', 'Input matrix D_k_km must be 3x3');
end

r_k = r_km + C_km' * d_k_km;
C_k = D_k_km * C_km;






