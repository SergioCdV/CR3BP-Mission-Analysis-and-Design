%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/05/21
% File: numerical_jacobian.m 
% Issue: 0 
% Validated: 

%% Numerical Hessian %%
% This file implements the function to compute the numerical Jacobian of
% any function using a 4 stencil central difference scheme

% Inputs: - handle function func, the function to derivate
%         - vector x, the point at which to evaluate the Hessian

% Outputs: - array H, the required Hessian

% Methods: 4th order central difference scheme

% New versions:

function [H] = numerical_hessian(func, x)
    % Constants 
    e = 1e-8;                               % Pertubation step 
    delta = eye(length(x));                 % Identity matrix / Kronecker delta
    stencil = 4;                            % Order of the stencil 
    K = [2 1 -1 -2];                        % Stencil coefficients
    
    % Preallocation 
    H = zeros(length(x), length(x));        % Hessian matrix
        
    % Main computation
    for i = 1:length(x) 
        for j = 1:length(x)
            if (i == j)
                % Perturbed state
                P = zeros(length(x),stencil);               % Preallocation of the perturbed vectors
                for k = 1:stencil
                    P(:,k) = x + K(k) * e * delta(:,i);     % Perturbed reference point
                end     
                
                % Hessian matrix
                H(i,i) = (-feval(func, P(:,4))-feval(func, P(:,1))      ...
                          +16*(feval(func, P(:,2))+feval(func, P(:,3))) ... 
                          -30*feval(func, x))/(12*e^2);                               % Hessian entry
            else
                % Perturbed state
                P(:,1) = x + e * delta(:,i) + e * delta(:,j);                         % Perturbed reference point
                P(:,2) = x + e * delta(:,i) - e * delta(:,j);                         % Perturbed reference point  
                P(:,3) = x - e * delta(:,i) + e * delta(:,j);                         % Perturbed reference point  
                P(:,4) = x - e * delta(:,i) - e * delta(:,j);                         % Perturbed reference point  
                
                % Hessian matrix
                H(i,j) = (feval(func, P(:,1))+feval(func, P(:,4))-feval(func, P(:,2))-feval(func, P(:,3)))/(4*e^2);                               
            end
        end
    end
end