%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(params, beta, t0, tf, s, u)
    % Minimum fuel
    M = 0; 
    L = sqrt(dot(u,u,1));

    % Minimum energy
%     M = 0; 
%     L = dot(u,u,1);

    % Minimum time
%     M = 0;
%     L = ones(1,size(s,2));
end