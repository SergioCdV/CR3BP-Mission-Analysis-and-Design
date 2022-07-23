%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/04/22
% File: CH_derivative.m 
% Issue: 0 
% Validated: 

%% Chebysev derivative %%
% This function allows to compute all Chebyshev polynomials derivatives of order n of both kinds,
% evaluated at the argument u. 

% Inputs: - string kind, specifying the kind of Chebyshev polynomials to
%           compute
%         - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated
%         - scalar degree, the degree of the derivative to be computed

% Outpus: - vector Pn, containing the evaluated Chebyshev polynomials
%           derivatives

function [B] = CH_derivative(kind, order, u, degree)
    % Switch the derivative order
    switch (degree)
        case 1
            B = dchebyshev(kind, order, u);
        case 2
            B = ddchebyshev(kind, order, u);
        otherwise
            error('A higher-order Chebyshev polynomial derivative is required, but has not been implemented')
    end
end

%% Auxiliary functions 
% First order basis of the Chebyshev tangent space
function [dPn] = dchebyshev(kind, order, u)
    % Preallocation of the polynomials and its derivatives
    Pn = CH_basis(kind,order,u);
    dPn = zeros(order+1,length(u)); 

    % Main computation 
    switch (kind)
        case 'first'
            dPn(1,:) = zeros(1,length(u));      % Initialization of the Chebyshev polynomials of the first kind
            dPn(2,:) = ones(1,length(u));       % Initialization of the Chebyshev polynomials of the first kind

        case 'second'
            dPn(1,:) = zeros(1,length(u));      % Initialization of the Chebyshev polynomials of the second kind
            dPn(2,:) = 2*ones(1,length(u));     % Initialization of the Chebyshev polynomials of the second kind

        otherwise
            error('No valid kind of polynomials was selected'); 
    end
  
    % Chebyshev polynomials derivatives 
    for i = 2:order
        dPn(i+1,:) = 2*Pn(i,:)+2*u.*dPn(i,:)-dPn(i-1,:);  
    end
end

% Second order basis of the Chebyshev tangent space
function [ddPn] = ddchebyshev(kind, order, u)
    % Preallocation of the polynomials
    dPn = dchebyshev(kind, order, u);
    ddPn = zeros(order+1,length(u));  

    % Main computation 
    switch (kind)
        case 'first'
            ddPn(1,:) = zeros(1,length(u));                    % Initialization of the Chebyshev polynomials of the first kind
            ddPn(2,:) = zeros(1,length(u));                    % Initialization of the Chebyshev polynomials of the first kind

        case 'second'
            ddPn(1,:) = zeros(1,length(u));                    % Initialization of the Chebyshev polynomials of the second kind
            ddPn(2,:) = zeros(1,length(u));                    % Initialization of the Chebyshev polynomials of the second kind

        otherwise
            error('No valid kind of polynomials was selected'); 
    end
  
    % Chebyshev polynomials derivatives 
    for i = 2:order
        ddPn(i+1,:) = 2*u.*ddPn(i,:)+4*dPn(i,:)-ddPn(i-1,:);  
    end
    
end