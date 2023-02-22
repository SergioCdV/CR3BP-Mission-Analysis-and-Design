%% Project:  
% Sergio Cuevas del Valle
% Date: 09/05/22
% File: OB_nodes
% Issue: 0 
% Validated: 

%% Orthogonal Bernstein nodes %%
% This scripts provides the function to compute the Legendre-Gauss nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the degree of the Legendre polynomial of interest

% Output: - vector y, the orthogonal Bernstein nodes of interest

function [y] = OB_nodes(N)
    % Polynomial orders
    N = N-1;
    N(2) = N(1)+1; 
    N(3) = N(1)+2;
    
    % Initial guess using Gauss-Radau nodes
    y = linspace(0, 1, N(2)).';
        
    % Compute the zeros of the N+1 Legendre Polynomial using the recursion relation and the Newton-Raphson method
    y0 = 2;         % Convergence value
    tol = 1e-10;    % Tolerance

    iter = 1;       % Initial iteration 
    maxIter = 1e4;  % Maximum number of iterations

    % Newton-Rhapson method
    while (max(abs(y-y0)) > tol && iter < maxIter)   
        L = OB_basis(N(2),y).';
        dL = OB_derivative(N(2),y,1).';

        dy = -L(:,N(3))./dL(:,N(3));
        y0 = y; 
        y = y0+dy;

        iter = iter+1;
    end
end