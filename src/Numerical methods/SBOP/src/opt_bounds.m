%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Boundary conditions %% 
% Function to compute the initial trajectory guess given the boundary
% conditions of the spacecraft in non-dimensional units

% Inputs: - class Problem, defining the problem of interest
%         - vector n, the order of the approximating polynomial function 
%         - scalar B, the size of the beta variables

% Outputs: - vector P_lb, vector of lower bounds
%          - vector P_ub, vector of upper bounds

function [P_lb, P_ub] = opt_bounds(Problem, n, B)
    % Constants 
    StateDim = Problem.StateDim; 

    % Upper and lower bounds
    [LB, UB] = Problem.Bounds(); 

    if (isempty(LB))
        P_lb = -Inf * ones(1,2 + (max(n)+1) + B).';
    else
        P_lb = [repmat(LB(1:StateDim), 1, max(n)+1) LB(StateDim+1:end)].';
    end

    if (isempty(UB))
        P_ub = Inf * ones(1,2 + (max(n)+1) + B).';
    else
        P_ub = [repmat(UB(1:StateDim), 1, max(n)+1) UB(StateDim+1:end)].';
    end
end