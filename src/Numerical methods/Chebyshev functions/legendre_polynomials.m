%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 02/07/21
% File: legendre_polynomials.m 
% Issue: 0 
% Validated: 

%% Legendre polynomials %%
% This function allows to compute all Legendre polynomials of order n of both kinds,
% evaluated at the argument x. 

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated

% Outpus: - vector Pn, containing the evaluated Legendre polynomials 

function [Pn] = legendre_polynomials(order, u)
    %Preallocation of the polynomials 
    Pn = zeros(order,1); 

    for i = 1:order
        n = i-1;
        for k = 0:n
            Pn(i) = Pn(i)+(1/2^n)*(factorial(n)/(factorial(k)*factorial(n-k)))^2*(u+1)^(n-k)*(u-1)^k; 
        end
    end
end
