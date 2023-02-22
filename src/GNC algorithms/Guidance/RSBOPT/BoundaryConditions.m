%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Boundary Conditions function %% 
% Function implementation of the boundary conditions definition

function [s0, sf] = BoundaryConditions(initial, final, beta, t0, tf)
    s0 = initial; 
    sf = final;
    sf(2) = beta(1);
end