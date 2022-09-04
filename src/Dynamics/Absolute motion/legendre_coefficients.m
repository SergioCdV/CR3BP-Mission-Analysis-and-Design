%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: legendre_coefficients.m 
% Issue: 0 
% Validated: 

%% Legendre coefficients %%
% For a given gravitational parameter mu and a colinear libration point, this function computes
% the associated Legendre polynomials coefficients

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - scalar libration point L, 1 for L1, 2 for L2 and 3 for L3
%         - scalar gamma, distance from the libration point to the least
%           massive primary
%         - scalar order, defining the order up to which the coefficients
%           will be computed

% Outputs: - the array c, containing the Legendre coefficients

% New versions: 

function [c] = legendre_coefficients(mu, L, gamma, order)
    % Preallocation 
    c = zeros(1, order);
    
    % Main computation 
    switch (L)
        case 1
            for i = 2:order
               c(i) = (1/gamma^3)*(mu+(-1)^i*((1-mu)*gamma^(i+1))/(1-gamma)^(i+1));
            end
        case 2
            for i = 2:order
               c(i) = ((-1)^i/gamma^3)*(mu+((1-mu)*gamma^(i+1))/(1+gamma)^(i+1));
            end
        case 3
            for i = 2:order
               c(i) = ((-1)^i/gamma^3)*(1-mu+((mu*gamma^(i+1))/(1+gamma)^(i+1)));
            end
        otherwise
            error('No valid Lagrange point was selected');
    end
end