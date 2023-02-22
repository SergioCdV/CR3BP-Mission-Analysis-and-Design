%% Project:  
% Sergio Cuevas del Valle
% Date: 06/11/22
% File: LG_weights
% Issue: 0 
% Validated: 

%% Legendre-Gauss weights %%
% This scripts provides the function to compute the Legendre-Gauss weights
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the order of the quadrature

% Output: - vector w, the Legendre weights of interest
%         - vector tau, the Legendre nodes of Pn

function [w, tau] = LG_weights(N)    
    beta = .5./sqrt(1-(2*(1:N)).^(-2));         % 3-term recurrence coeffs
    T = diag(beta,1) + diag(beta,-1);           % Jacobi matrix
    [V,D] = eig(T);                             % Eigenvalue decomposition
    x = diag(D);                                % Eigenvalue matrix  
    [tau,i] = sort(x);                          % Legendre nodes
    w = 2*V(1,i).^2;                            % Weights
end