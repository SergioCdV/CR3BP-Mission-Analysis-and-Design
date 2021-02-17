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

% Output: - the index s, containing information about close bifurcations
%           around the solution associated with the STM.

%New versions: use symplecticity to correct computational error.

function [s, state] = henon_stability(STM)
    %Compute STM eigenspectrum 
    [~, eig, flag] = eigs(STM);
            
    %Compute Henon Stability Index
    uLambda = eig(1);               %Greatest eigenvalue, corresponding to the unstable direction of the manifold
    sLambda = 1/uLambda;            %Reciprocal of the unstable eigenvalue, corresponding to the stable direction
    s = (1/2)*(sLambda+uLambda);    %Henon Stability Index
    
    %Sanity check 
    if (flag == 1)
        state = false;              %Eigen-decomposition of the STM failed
    else
        state = true;               %Eigen-decomposition of the SMT succeed
    end
end
