%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/12/24
% File: legendre_coefficients.m 
% Issue: 0 
% Validated: 

%% Legendre coefficients %%
% For a given gravitational parameter mu and a colinear libration point, this function computes
% the associated Legendre polynomials coefficients

% Inputs: - double mu, the reduced gravitational parameter of the system
%         - double libration point L, 1 for L1, 2 for L2 and 3 for L3
%         - double gamma, distance from the libration point to the least
%           massive primary
%         - double order, defining the order up to which the coefficients
%           will be computed

% Outputs: - array c [1 x order+1], containing the Legendre coefficients

% New versions: 

function [c] = LegendreCoefficients(mu, L, gamma, order)
    % Preallocation 
    c = zeros(1, order+1);
    order_array = 2:order;
    
    % Main computation 
    switch (L)
        case 1
            c(3:end) = ( mu + (-1).^order_array .* ( (1 - mu) * gamma.^(order_array + 1) ) ./ (1 - gamma).^(order_array + 1) ) / gamma^3;

        case 2
            c(3:end) = ( (-1).^order_array / gamma^3 ) .* ( mu + ( (1 - mu) * gamma.^(order_array + 1) ) ./ (1 + gamma).^(order_array + 1) );

        case 3
            c(3:end) = ( (-1).^order_array / gamma^3 ) .* ( 1 - mu + ( mu * gamma.^(order_array + 1) ./ (1 + gamma).^(order_array + 1) ) );

        otherwise
            error('No valid Lagrange point was selected');
    end
end