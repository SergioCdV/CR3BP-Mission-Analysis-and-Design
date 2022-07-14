%% Project: 
% Date: 30/01/2022

%% Bernstein derivative %% 
% Function to compute the basis of the tangent space of the Bernstein polynomials 

% Inputs:  - scalar n, the cardinal of the polynomial basis
%          - vector tau, the parametrization of the basis
%          - scalar order, the order of the derivative required

% Outputs: - array B, the basis of the required tangent space, of size n+1 x length(tau)

function [B] = bernstein_derivative(n, tau, order)
    % Switch the derivative order
    switch (order)
        case 1
            B = dbernstein(n, tau);
        case 2
            B = ddbernstein(n, tau);
        otherwise
            error('A higher-order Bernstein polynomial derivative is required, but has not been implemented')
    end
end

%% Auxiliary functions 
% First order basis of the Bernstein tangent space
function [dB] = dbernstein(n, tau)
    % Preallocation for speed 
    dB = zeros(n+1,length(tau));

    if (n ~= 0)
        B = bernstein_basis(n-1,tau);
    
        % Computation of the derivative basis 
        for i = 0:n
            if (i == 0)
                dB(i+1,:) = -n*B(i+1,:);
            elseif (i == n)
                dB(i+1,:) = n*B(i,:);
            else
                dB(i+1,:) = n*(B(i,:)-B(i+1,:));
            end
        end
    end
end

% Second order basis of the Bernstein tangent space 
function [dB] = ddbernstein(n, tau)
    % Preallocation for speed 
    dB = zeros(n+1,length(tau));

    if (n ~= 0)
        B = dbernstein(n-1,tau);
    
        % Computation of the derivative basis 
        for i = 0:n
            if (i == 0)
                dB(i+1,:) = -n*B(i+1,:);
            elseif (i == n)
                dB(i+1,:) = n*B(i,:);
            else
                dB(i+1,:) = n*(B(i,:)-B(i+1,:));
            end
        end
    end
end