function [A] = xMat(a)

% function xMat
% Keith Leung 2013
%
% Turns a 3d vector into a cross-product matrix such that
% a x b = Ab
% 
% Input: a - 3d vector
% output: A - cross-product matrix

if(nargin ~= 1)
    error('xMat:nargChk', 'xMat takes 1 input only');
end
if(length(a(:,1)) ~= 3 || length(a(1,:)) ~= 1)
    error('xMat:sizeChk', 'Input vector must be 3x1');
end

A = [0 -a(3) a(2);
     a(3) 0 -a(1);
     -a(2) a(1) 0];