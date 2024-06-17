%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 17/06/24
% File: test_1.m 
% Issue: 0 
% Validated: 

%% Test 1 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 1 is concerned with the generation of the basic objects of the
% library

%% Input parameters 
% System characteristics
mu = 0.0121505856;                       % Mass parameter for the Earth-Moon system

% Create a CR3BP system 
EarthMoon = src.Systems.CR3BPSystem( 'mu', mu );