%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/09/22
% File: HT_derivative.m 
% Issue: 0 
% Validated:

%% Hermite polynomials derivative %%
% This function allows to compute all Hermite polynomials derivatives of order n,
% evaluated at the argument u

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated
%         - scalar degree, the degree of the derivative to be computed

% Outpus: - vector Pn, containing the evaluated Hermite polynomials
%           derivatives

function [B] = HT_derivative(order, u, degree)
    % Switch the derivative order
    switch (degree)
        case 1
            B = dhermite(order, u);
        case 2
            B = ddhermite(order, u);
        otherwise
            error('A higher-order Hermite polynomial derivative is required, but has not been implemented')
    end
end

%% Auxiliary functions 
% First order basis of the Hermite tangent space
function [dPn] = dhermite(order, u)
    % Preallocation of the polynomials and its derivatives
    Pn = HT_basis(order,u);
    dPn = zeros(order+1,length(u)); 

    % Bonnet's formula 
    for i = 2:order+1
        n = i-1;
        dPn(i,:) = 2*n*Pn(i-1,:); 
    end
end

% Second order basis of the Hermite tangent space
function [ddPn] = ddhermite(order, u)
    % Preallocation of the polynomials and its derivatives
    Pn = HT_basis(order,u);
    dPn = HT_derivative(order,u,1); 
    ddPn = zeros(order+1,length(u)); 

    % Bonnet's formula 
    for i = 1:order+1
        n = i-1;
        ddPn(i,:) = 2*u.*dPn(i,:)-2*n*Pn(i,:); 
    end
end