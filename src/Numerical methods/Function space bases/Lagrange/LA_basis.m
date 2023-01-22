%% Project: 
% Date: 08/11/22

%% Lagrange basis %%
% This function allows to compute all Lagrange polynomials of order n,
% evaluated at the argument u. 

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated

% Outpus: - vector Pn, containing the evaluated Lagrange polynomials 

function [Pn] = LA_basis(order, u)
    % Preallocation of the polynomials 
    Pn = zeros(order+1,length(u)); 

    % Definition of Lagrange polynomials
    for i = 1:order+1
        Pn(i,:) = ones(1,length(u));
        for j = 1:length(u)
            if (i ~= j)
                Pn(i,:) = Pn(i,:).*(u-u(j))/(u(i)-u(j));
            end
        end
    end
end
