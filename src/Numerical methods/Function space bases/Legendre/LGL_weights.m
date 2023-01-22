%% Project:  
% Sergio Cuevas del Valle
% Date: 06/11/22
% File: LG_weights
% Issue: 0 
% Validated: 

%% Legendre-Gauss-Lobatto weights %%
% This scripts provides the function to compute the Legendre-Gauss-Lobatto weights
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the order of the quadrature

% Output: - vector w, the Legendre weights of interest
%         - vector tau, the Legendre nodes of Pn
%         - matrix D, the differentiation matrix

function [W, tau, D] = LGL_weights(N)    
    % Truncation + 1
    N1 = N+1;

    % CGL nodes
    xc = cos(pi*(0:N)/N)';

    % Uniform nodes
    xu = linspace(-1,1,N1).';

    % Make a close first guess to reduce iterations
    if (N < 3)
        tau = xc;
    else
        tau = xc + sin(pi*xu)./(4*N);
    end

    % Preallocation of the Legendre Vandermonde Matrix
    P = zeros(N1);

    % Compute P_(N) using the recursion relation. Compute its first and second derivatives and  update x using the Newton-Raphson method.
    xold = 2;

    while (max(abs(tau-xold)) > eps)
        % Initialization
        xold = tau;
        P(:,1) = 1;    
        P(:,2) = tau;
        
        for k = 2:N
            P(:,k+1)=( (2*k-1)*tau.*P(:,k)-(k-1)*P(:,k-1) )/k;
        end
         
        tau = xold-( tau.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
    end

    % Legendre basis 
    tau = flip(tau).';
    L = LG_basis(N, tau);

    % Quadrature weights
    W = 2./(N*N1*L(end,:).^2);

    % Differentiation matrix
    D = zeros(N+1);
    for i = 0:N 
        for j = 0:N
            if (i ~= j)
                D(i+1,j+1) = (L(end,i+1)/L(end,j+1))/(tau(i+1)-tau(j+1));
            elseif (i == j && j == N)
                D(i+1,j+1) = (N*(N+1)/4);
            elseif (i == j && j == 0)
                D(i+1,j+1) = -(N*(N+1)/4);
            else 
                D(i+1,j+1) = 0; 
            end
        end
    end

    tau = tau.';
end