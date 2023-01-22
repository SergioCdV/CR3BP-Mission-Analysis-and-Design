%% Project: 
% Date: 30/04/22

%% Legendre derivative %%
% This function allows to compute all Legendre polynomials derivatives of order n,
% evaluated at the argument u. 

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated
%         - scalar degree, the degree of the derivative to be computed

% Outpus: - vector Pn, containing the evaluated Legendre polynomials
%           derivatives

function [B] = LG_derivative(order, u, degree)
    % Switch the derivative order
    switch (degree)
        case 1
            B = dlegendre(order, u);
        case 2
            B = ddlegendre(order, u);
        otherwise
            error('A higher-order Legendre polynomial derivative is required, but has not been implemented')
    end
end

%% Auxiliary functions 
% First order basis of the Legendre tangent space
function [dPn] = dlegendre(order, u)
    % Preallocation of the polynomials and its derivatives
    Pn = LG_basis(order,u);
    dPn = zeros(order+1,length(u)); 

    % Initialization of the polynomials 
    dPn(1,:) = zeros(1,length(u)); 
    dPn(2,:) = ones(1,length(u)); 

    % Bonnet's formula 
    for i = 3:order+1
        n = i-1;
        dPn(i,:) = ((2*n-1)*(Pn(i-1,:)+u.*dPn(i-1,:))-(n-1)*dPn(i-2,:))/n; 
    end
end

% Second order basis of the Legendre tangent space
function [ddPn] = ddlegendre(order, u)
    % Preallocation of the polynomials and its derivatives
    dPn = LG_derivative(order,u,1); 
    ddPn = zeros(order+1,length(u));

    % Initialization of the polynomials 
    ddPn(1,:) = zeros(1,length(u)); 
    ddPn(2,:) = zeros(1,length(u));  

    % Bonnet's formula 
    for i = 3:order+1
        n = i-1;
        ddPn(i,:) = ((2*n-1)*(2*dPn(i-1,:)+u.*ddPn(i-1,:))-(n-1)*ddPn(i-2,:))/n; 
    end
end