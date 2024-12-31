%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 19/12/24
% File: AugmentedPotential.m 
% Issue: 0 
% Validated: 

%% Augmented potential %%
% For a given gravitational parameter mu and position vectors r, this function computes the 
% augmented potential function associated with that input position vector

% Inputs: - double mu, the reduced gravitational parameter of the system 
%         - array r [12xN], containing the synodic co-orbital state vectors

% Outputs: - vector U [N], the augmented potential function associated of
%            the co-orbital problem

% New versions:

function [U] = AugmentedPotential(mu, r)
    % Obtain synodic coordinates 
    x = r(1,:);                         % Synodic x coordinate
    y = r(2,:);                         % Synodic y coordinate 
    
    % Augmented potential function
    U = src.Systems.CR3BPSystem.CoPotentialFunction(mu, r);
    U = -0.5 * (x.^2 + y.^2) + U;
end