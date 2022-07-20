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
    % Initial guess 
    y = (1-(3*(N-2))/(8*(N-1)^3))*cos(pi*(4*(2:N-1)-3)/(4*(N-1)+1)).';

    % Polynomial orders
    N = N-1;
    
    % Preallocation of the Legendre polynomials
    L = zeros(length(y),N+1);
    dL = L; 
    
    % Compute the zeros of the N+1 Legendre Polynomial using the recursion relation and the Newton-Raphson method
    y0 = 2;         % Convergence value

    iter = 1;       % Initial iteration
    maxIter = 1e4;  % Maximum number of iterations

    % Newton-Rhapson method
    while (max(abs(y-y0)) > eps & iter < maxIter)   
        for i = 1:length(y)
            L(i,:) = LG_derivative(N,y(i),1).';
            dL(i,:) = LG_derivative(N,y(i),2).';
        end

        dy = -L(:,end)./dL(:,end);
        y0 = y; 
        y = y0+dy;

        iter = iter+1; 
    end
    
    % Linear map from [-1,1] to [a,b]
    y = flip(y);
    y = [-1; y; 1];
end