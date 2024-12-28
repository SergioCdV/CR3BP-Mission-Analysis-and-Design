%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/12/24
% File: VariationalEquationsCR3BP.m 
% Issue: 0 
% Validated: 

%% Linear variational Equations of the co-orbital CR3BP dynamics %%
% This function contains the first order variational equations of the co-orbital CR3BP

% Input: -
% Output: - array ds [36xN], the differential equations of the problem 

function [ds] = VariationalEquationsCoCR3BP(t, j, s, params)
    % State transition matrix of the system 
    idx = size(s,1) - params(1)^2 + 1;

    Phi = s(idx:end,:);                           
    Phi = reshape(Phi, params(1), params(1) * size(s,2));

    % Reference trajectory 
    x_ref = s(1:idx, :);                          % Reference trajectory

    % Compute the first order variational equations 
    J = src.Systems.CoCR3BPSystem.JacobianCoCR3BP(params(2), x_ref);

    % Differential system 
    ds = J * Phi; 
    ds = reshape(ds, [], 1);
end