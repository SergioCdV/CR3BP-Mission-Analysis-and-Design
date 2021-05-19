%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 19/05/21
% File: OptControl_guidance.m 
% Issue: 0 
% Validated: 19/05/21

%% Optimal control guidance %%
% This script contains the function to compute the guidance law by means of a NPL formulation of the multi-impulsive
% scheme.

% Inputs: - string dynamics, to select the dsired type of APF algorithm
%         - boolean safe_corridor, to activate/desactivate the safety APF
%         - vector state, the state of the object from which the guidance
%           law is computed
%         - array obstacles_state, indicating the position of the obstacles
%           (if any) to avoid
%         - scalar dt, the time step taken in the simulation 
%         - vector Sg0, the initial conditions of the desired guidance law

% Output: - the index s, containing information about close bifurcations around the solution associated with the STM.

% New versions: 

function [] = optimal_control()
end