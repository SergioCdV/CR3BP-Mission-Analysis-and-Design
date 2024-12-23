%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/12/24
% File: test_2.m 
% Issue: 0 
% Validated: 

%% Test 2 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 2 is concerned with the generation of the basic objects representing
% periodic orbits in a CR3BP system

%% Basic system 
% Create a CR3BP system 
EarthMoon = src.Systems.EarthMoon();

r = EarthMoon.ZeroVelocitySurface(EarthMoon.mu, 3.1776, false);

%% Evaluation of the Jacobian of the absolute dynamics
s = [r; zeros(3, size(r,2))];
J = EarthMoon.JacobianCR3BP(EarthMoon.mu, s);

%% Evaluation of the vector field of the absolute dynamics
s = [r; zeros(3, size(r,2))];
ds = EarthMoon.Dynamics(0, 0, s, zeros(3,size(s,2)), EarthMoon.params);
