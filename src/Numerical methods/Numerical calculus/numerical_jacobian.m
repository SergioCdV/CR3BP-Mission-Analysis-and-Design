%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/05/21
% File: numerical_jacobian.m 
% Issue: 0 
% Validated: 

%% Numerical Jacobian %%
% This file implements the function to compute the numerical Jacobian of
% any function using a 4 stencil central difference scheme

% Inputs: - handle function func, the function to derivate.
%         - vector x, the point at which to evaluate the Jacobian. 

% Outputs: - array J, the required Jacobian.

% Methods: 4th order central difference scheme.

% New versions:

function [J] = numerical_jacobian(dim, func, x)
    %Constants 
    e = 1e-6;                               %Pertubation step 
    delta = eye(length(x));                 %Identity matrix / Kronecker delta
    stencil = 4;                            %Order of the stencil 
    K = [2 1 -1 -2];                        %Stencil coefficients
    
    %Preallocation 
    J = zeros(dim, length(x));              %Preallocation of the Jacobian 
    
    %Main computation
    for i = 1:length(x) 
        P = zeros(length(x),stencil);       %Preallocation of the perturbed vectors
        for j = 1:stencil
            P(:,j) = x + K(j) * e * delta(:,i);
        end
        J(:,i) = (feval(func, P(:,4))-feval(func, P(:,1))+8*(feval(func, P(:,2))-feval(func, P(:,3))))/(12*e);             
    end
end