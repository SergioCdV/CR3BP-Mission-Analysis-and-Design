%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: henon_stability.m 
% Issue: 0 
% Validated: 10/12/20

%% Henon Stability Index %%
% This script contains the function to compute the Henon stability index
% associated with a dynamical solution of the problem.

% Inputs: - state transition matrix STM, whose eigenvalues are to be analyzed.

% Output: - the index s, containing information about close bifurcations around the solution associated with the STM.

%New versions: use symplecticity to correct computational error.

function [s, state] = henon_stability(STM)
    %Compute STM eigenspectrum 
    [~, eig, flag] = eigs(STM);
               
    %Henon stability indexes
    s(1) = (1/2)*(eig(1,1)+eig(6,6));    %Sum of the reciprocal pair
    s(2) = (1/2)*(eig(2,2)+eig(3,3));    %Sum of the neutrally stable pair
    s(3) = (1/2)*(eig(4,4)+eig(5,5));    %Sum of the remaining pair
    
    %Sanity check 
    if (flag == 1)
        state = false;      %Eigen-decomposition of the STM failed
    else
        state = true;       %Eigen-decomposition of the STM succeed
    end
end
