function [d, rot] = velocities2Displacement(v, w, dt)

% Function velocities2Displacement
% Keith Leung 2013
%
% Determines displacement from constant velocities over dt
%
% Inputs:
% v - translational velocity
% w - rotational velocity
% dt - timestep
%
% Outputs:
% d - position displacement
% rot - rotation

d = v * dt;
rot = w * dt;
