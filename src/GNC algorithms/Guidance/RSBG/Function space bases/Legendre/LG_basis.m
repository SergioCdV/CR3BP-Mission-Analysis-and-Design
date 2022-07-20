%% Project: 
% Date: 30/04/22

%% Legendre basis %%
% This function allows to compute all Legendre polynomials of order n,
% evaluated at the argument u. 

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated

% Outpus: - vector Pn, containing the evaluated Legendre polynomials 

function [Pn] = LG_basis(order, u)
    % Preallocation of the polynomials 
    Pn = zeros(order+1,length(u)); 

    % Initialization of the polynomials 
    Pn(1,:) = ones(1,length(u)); 
    Pn(2,:) = u; 

    % Bonnet's formula 
    for i = 3:order+1
        n = i-1;
        Pn(i,:) = ((2*n-1)*u.*Pn(i-1,:)-(n-1)*Pn(i-2,:))/n; 
    end
end
