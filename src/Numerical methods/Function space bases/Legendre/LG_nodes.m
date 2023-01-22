%% Project:  
% Sergio Cuevas del Valle
% Date: 20/01/20
% File: LG_nodes
% Issue: 0 
% Validated: 

%% Legendre-Gauss nodes %%
% This scripts provides the function to compute the Legendre-Gauss nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the degree of the Legendre polynomial of interest 

% Output: - vector y, the Legendre nodes of interest

function [y] = LG_nodes(N)    
    beta = .5./sqrt(1-(2*(1:N)).^(-2));         % 3-term recurrence coeffs
    T = diag(beta,1) + diag(beta,-1);           % Jacobi matrix
    [~,D] = eig(T);                             % Eigenvalue decomposition
    x = diag(D);                                % Eigenvalue matrix  
    y = sort(x);                                % Legendre nodes
end