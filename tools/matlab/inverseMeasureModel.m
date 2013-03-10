function [p_m_i] = inverseMeasureModel(p_k_i, C_k_i, p_m_k)

% function inverseMeasureModel
% Keith Leung 2013
%
% Inverse measurement model
% 
% Inputs:
% p_k_i - vehicle position relative to inertial frame
% C_k_i - vehicle orientation; rotation from inertial to vehicle frame
% p_m_k - feature position relative to vehicle frame (i.e., measurement)
%
% Output:
% p_m_i - feature position in inertial frame
%

if(nargin ~= 3)
    error('inverseMeasureModel:nargChk', 'measureModel takes 1 input only');
end
if(length(p_k_i(:,1)) ~= 3 || length(p_k_i(1,:)) ~= 1)
    error('inverseMeasureModel:sizeChk', 'Input vector p_k_i must be 3x1');
end
if(length(C_k_i(:,1)) ~= 3 || length(C_k_i(1,:)) ~= 3)
    error('inverseMeasureModel:sizeChk', 'Input matrix C_k_i must be 3x3');
end
if(length(p_m_k(:,1)) ~= 3 || length(p_m_k(1,:)) ~= 1)
    error('inverseMeasureModel:sizeChk', 'Input vector p_m_k must be 3x1');
end

p_m_i = C_k_i'*p_m_k + p_k_i;
