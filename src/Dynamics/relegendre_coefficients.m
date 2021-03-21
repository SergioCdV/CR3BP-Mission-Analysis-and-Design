%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 21/03/20
% File: relegendre_coefficients.m 
% Issue: 0 
% Validated: 

%% Relative Legendre coefficients %%
% This function computes the Legendre coefficients of the Legendre expasion of the relative motion. 

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - vector r_t, the target spacecraft synodic position vector.
%         - scalar order, the order of the expansion

% Outputs: - vector cn, containing the Legendre coefficients up to order.

% Methods: . 

% New versions: .

function [cn] = relegendre_coefficients(mu, r_t, order)
    %Characteristics of the system 
    mu1 = 1-mu;                         %Reduced gravitational parameter of the first primary
    mu2 = mu;                           %Reduced gravitational parameter of the second primary
    R1 = [-mu; 0; 0];                   %Synodic position of the first primary
    R2 = [1-mu; 0; 0];                  %Synodic position of the second primary
    
    %Preallocation of the coefficients 
    cn = zeros(1, order);
    
    %Main computation
    for i = 2:order
        cn(i) = mu1/norm(R1-r_t)^(i+1)+mu2/norm(R2-r_t)^(i+1);
    end
end