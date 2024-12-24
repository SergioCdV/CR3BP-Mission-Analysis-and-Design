%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/12/24
% File: VariationalEquationsCR3BP.m 
% Issue: 0 
% Validated: 

%% Linear variational Equations of the CR3BP dynamics %%
% This function contains the first order variational equations of the CR3BP

% Input: -
% Output: - array ds [36xN], the differential equations of the problem 

function [ds] = VariationalEquationsCR3BP(t, j, s, params)
    % Arrange the STM
    s = reshape(s, params(1) + params(1)^2, []);        % Arrange the state
    x_ref = s(1:params(1), :);                          % Reference trajectory

    % State transition matrix of the system 
    Phi = s(params(1)+1:end,:);                           
    Phi = reshape(Phi, params(1), params(1) * size(x_ref,2));

    % Compute the first order variational equations 
    J = src.Systems.CR3BPSystem.JacobianCR3BP(params(2), x_ref);

    % Differential system 
    ds = J * Phi; 
    ds = reshape(ds, [], 1);
end