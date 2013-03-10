function [p_m_k] = measureModel(p_k_i, C_k_i, p_m_i)

% 3d point measurement model
% Keith Leung 2013
%
% Given sensor pose and a feature position, returns the relative 3d
% position measurement
%
% Inputs: 
% p_k_i - sensor 3d position relative to inertial frame
% C_k_i - sensor orientation relative to inertial frame expressed as a
%       rotation matrix
% p_m_i - feature 3d position relative to inertial frame
%
% Outputs:
% p_m_k - relative 3d position measurement in vehicle frame

if(nargin ~= 3)
    error('measureModel:nargChk', 'measureModel takes 1 input only');
end
if(length(p_k_i(:,1)) ~= 3 || length(p_k_i(1,:)) ~= 1)
    error('measureModel:sizeChk', 'Input vector p_k_i must be 3x1');
end
if(length(C_k_i(:,1)) ~= 3 || length(C_k_i(1,:)) ~= 3)
    error('measureModel:sizeChk', 'Input matrix C_k_i must be 3x3');
end
if(length(p_m_i(:,1)) ~= 3 || length(p_m_i(1,:)) ~= 1)
    error('measureModel:sizeChk', 'Input vector p_m_i must be 3x1');
end

p_m_k = C_k_i * (p_m_i - p_k_i);  

