function [D] = aa2Rmat(w)

% function euler2Rmat
% Keith Leung 2013
%
% Express an axis-angle rotation as a rotation matrix
%
% Input: w - axis-angle vector
% Output: D - rotation matrix

if(nargin ~= 1)
    error('aa2Rmat:nargChk', 'aa2Rmat takes 1 input only');
end
if(length(w(:,1)) ~= 3 || length(w(1,:)) ~= 1)
    error('aa2Rmat:sizeChk', 'Input vector must be 3x1');
end

mag = norm(w);
c = cos(mag);
s = sin(mag);
e = w/mag;
D = eye(3)*c + (1-c)*e*e' - s*xMat(e);
