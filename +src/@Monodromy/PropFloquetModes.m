%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/01/25
% File: PropFloquetModes.m 
% Issue: 0 

%% Propagate Floquet modes %%
% This script propagates the Floquet modes to a given epoch

% Inputs: - array V [nxn], the Floquet multipliers of the monodromy
%         - array E0 [nxn], the Floquet vectors of the matrix at the period
%         - array Phi [nxn], the STM at epoch
%         - scalar t, the epoch at which to propagate the Floquet modes
%         - scalar T, the period of the monodromy matrix 

% Outputs: - array E [nxn], the Floquet vectors of the matrix at epoch

function [E] = PropFloquetModes(V, E0, Phi, t, T)
    % Floquet multipliers matrix 
    J = diag( exp( -t / T * diag(V) ) );

    % Propagation 
    E = Phi * E0 * J;
end