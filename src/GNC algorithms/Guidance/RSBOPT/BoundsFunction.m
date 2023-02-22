%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Uppe an lower bounds function %% 
% Function implementation the definition of the upper na dlowe bounds for
% the problem

function [LB, UB] = BoundsFunction()
    % Upper and lower bounds for the problem first order state vector, initial time, final time and parameters
    LB = [-Inf -Inf -Inf 0 0 0];
    UB = [Inf Inf Inf 0.1 Inf 30];
end