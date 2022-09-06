%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 21/05/21
% File: chebyshev.m 
% Issue: 0 
% Validated: 

%% Chebysev polynomials %%
% This function allows to compute all Chebyshev polynomials of order n of both kinds,
% evaluated at the argument u

% Inputs: - string kind, specifying the kind of Chebyshev polynomials to
%           compute
%         - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated

% Outpus: - vector Pn, containing the evaluated Chebyshev polynomials 

function [Pn] = chebyshev(kind, order, u)
    % Preallocation of the polynomials 
    Pn = zeros(order,1); 

    % Main computation 
    switch (kind)
        case 'first'
            Pn(1) = 1;                    % Initialization of the Chebyshev polynomials of the first kind
            Pn(2) = u;                    % Initialization of the Chebyshev polynomials of the first kind
        case 'second'
            Pn(1) = 1;                    % Initialization of the Chebyshev polynomials of the second kind
            Pn(2) = 2*u;                  % Initialization of the Chebyshev polynomials of the second kind
        otherwise
            error('No valid kind of polynomials was selected'); 
    end
  
    for i = 2:order-1
        Pn(i+1) = 2*u*Pn(i)-Pn(i-1);      % Chebyshev polynomials
    end
end
