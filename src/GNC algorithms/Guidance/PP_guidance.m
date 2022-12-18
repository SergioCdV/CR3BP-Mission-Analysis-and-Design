%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: setup_path.m 
% Issue: 0 
% Validated: 

%% Prussing's Pruner %%
% This scripts provides the function to prune an impulsive sequence based
% on primer vector theory for linear systems

function [dV] = PP_guidance(dV)
    % Compute the norm of the impulses 
    V = sqrt(dot(dV,dV,1)).'; 

    % Compute the alpha coefficients
    A = dV;                               % Linear coefficients  
    b = -rand.*dV(:,1);                   % Random constraint on alpha
    alpha = A\b;

    % Proceed with the pruner
    if (norm(alpha) == 0)
        % Not reduceable 
        dV = dV;
    else
        % Compute the beta coefficients 
        beta = alpha./V; 
    
        % Compute the mu coefficients 
        mu = V-alpha/max(beta);
    
        % Compute the new reduced sequence 
        dV = (mu./V).'.*dV;
    end
end