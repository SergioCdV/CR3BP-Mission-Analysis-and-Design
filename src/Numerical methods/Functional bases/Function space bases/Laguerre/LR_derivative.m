%% Project: 
% Date: 30/04/22

%% Laguerre derivative %%
% This function allows to compute all Laguerre polynomials derivatives of order n,
% evaluated at the argument u. 

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated
%         - scalar degree, the degree of the derivative to be computed

% Outpus: - vector Pn, containing the evaluated Laguerre polynomials
%           derivatives

function [B] = LR_derivative(order, u, degree)
    % Switch the derivative order
    switch (degree)
        case 1
            B = dlaguerre(order, u);
        case 2
            B = ddlaguerre(order, u);
        otherwise
            error('A higher-order Laguerre polynomial derivative is required, but has not been implemented')
    end
end

%% Auxiliary functions 
% First order basis of the Laguerre tangent space
function [dPn] = dlaguerre(order, u)
    % Preallocation of the polynomials and its derivatives
    Pn = LR_basis(order,u);
    dPn = zeros(order+1,length(u)); 

    % Initialization of the polynomials 
    dPn(1,:) = zeros(1,length(u)); 
    dPn(2,:) = -ones(1,length(u)); 

    % Bonnet's formula 
    for i = 2:order
        n = i-1;
        dPn(i+1,:) = ((2*n-1-u).*dPn(i,:)-Pn(i,:)-n*dPn(i-1,:))/(n+1); 
    end
end

% Second order basis of the Legendre tangent space
function [ddPn] = ddlaguerre(order, u)
    % Preallocation of the polynomials and its derivatives
    dPn = LR_derivative(order,u,1); 
    ddPn = zeros(order+1,length(u));

    % Initialization of the polynomials 
    ddPn(1,:) = zeros(1,length(u)); 
    ddPn(2,:) = zeros(1,length(u));  

    % Bonnet's formula 
    for i = 2:order
        n = i-1;
        ddPn(i+1,:) = ((2*n-1-u).*ddPn(i,:)-2*dPn(i,:)-n*ddPn(i-1,:))/(n+1); 
    end
end