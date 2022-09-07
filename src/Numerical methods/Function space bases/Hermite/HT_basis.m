%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/09/22
% File: HT_basis.m 
% Issue: 0 
% Validated:

%% Hermite basis %%
% This function allows to compute all physicist Hermite polynomials of order n,
% evaluated at the argument u

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated

% Outpus: - vector Pn, containing the evaluated Hermite polynomials 

function [Pn] = HT_basis(order, u)
    % Preallocation of the polynomials 
    Pn = zeros(order+1,length(u)); 

    % Initialization of the polynomials 
    Pn(1,:) = ones(1,length(u));
    Pn(2,:) = 2*u;

    % Bonnet's formula 
    for i = 2:order
        n = i-1;
        Pn(i+1,:) = 2*u.*Pn(i,:)-2*n*Pn(i-1,:); 
    end
end
