function [Sp_m_i] = inverseMeasureUncertainty(p_k_i, C_k_i, Sp_m_k)

% function inverseMeasureModel
% Keith Leung 2013
%
% Inverse measurement model
% 
% Inputs:
% p_k_i - vehicle position relative to inertial frame
% C_k_i - vehicle orientation; rotation from inertial to vehicle frame
% Sp_m_k - feature position uncertiany in vehicle frame
%
% Output:
% Sp_m_i - feature position uncertainty 
%

if(nargin ~= 3)
    error('inverseMeasureUncertainty:nargChk', 'measureModel takes 1 input only');
end
if(length(p_k_i(:,1)) ~= 3 || length(p_k_i(1,:)) ~= 1)
    error('inverseMeasureUncertainty:sizeChk', 'Input vector p_k_i must be 3x1');
end
if(length(C_k_i(:,1)) ~= 3 || length(C_k_i(1,:)) ~= 3)
    error('inverseMeasureUncertainty:sizeChk', 'Input matrix C_k_i must be 3x3');
end
if(length(Sp_m_k(:,1)) ~= 3 || length(Sp_m_k(1,:)) ~= 3)
    error('inverseMeasureUncertainty:sizeChk', 'Input covariance matrix Sp_m_k must be 3x3');
end

Sp_m_i = C_k_i' * Sp_m_k * C_k_i;
