%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/01/25
% File: FloquetAnalysis.m 
% Issue: 0 

%% Floquet Analysis %%
% This script computes the Floquet multipliers and vectors of a given
% Monodromy matrix

% Inputs: - array lambda [nxnN], the eigenvalues of the matrix
%         - array v [nxnN], the eigenvectors of the matrix
%         - scalar T, the period of the monodromy matrix 

% Outputs: - array V [nxnN], the Floquet multipliers of the matrix
%          - array E [nxnN], the Floquet vectors of the matrix

function [V, E] = FloquetAnalysis(lambda, v, T)
    % Floquet exponents
    V = log(lambda) / T;

    % Floquet modes
    E = v;                  % The initial conditions of the Floquet modes
end