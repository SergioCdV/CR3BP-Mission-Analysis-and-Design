%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: cr3bp_vecEquations.m 
% Issue: 0 
% Validated: 

%% CR3BP Vectorized Dynamics %%
% This function contains the vectorized dynamical model of the problem, 
% useful to integrate large initial conditions batches

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - scalar direction (in binary format, 1 or -1), indicating the
%           time integration direction: 1 for forward integration, -1 for
%           backward integration
%         - boolean flagVar, true for dyanmics and STM integration, 
%           false for only dynamical integration
%         - scalar t, a reference epoch
%         - vector s, containing in an Nx1 array the phase space vector,
%           possibly augmented with the state transition matrix at time t

% Outputs: - vector dr, the differential vector field, which will include
%            the phase space trajectory and the STM integrated state when
%            flagVar is true

% Methods: non-dimensional CR3BP dynamics in the synodic frame

% New versions: include STM vectorized integration

function [dr] = cr3bp_vecEquations(mu, direction, flagVar, t, s)    
    % Define the initial phase space vector
    x = s(1,:,:,:);         % Synodyc x coordinate
    y = s(2,:,:,:);         % Synodyc y coordinate 
    z = s(3,:,:,:);         % Synodyc z coordinate 
    V = s(4:6,:,:,:);       % Synodic velocity vector
    
    % Relevant system parameters
    mu1 = (1-mu)*ones(size(V,2),size(V,2),size(V,2));   % First primary normalized position
    mu2 = mu*ones(size(V,2),size(V,2),size(V,2));       % Second primary normalized position
    r1 = [(x+mu2); y; z];                               % Relative position vector to the first primary
    r2 = [(x-mu1); y; z];                               % Relative position vector to the secondary primary
    R1 = norm(r1);                                      % Relative distance to the first primary
    R2 = norm(r2);                                      % Relative distance to the second primary
    
    %Compute the time flow of the system (depending on the time direction)
    if (direction == -1)
        %Inertial acceleration
        inAcc = [x-2*V(2,:,:,:); y+2*V(1,:,:,:); zeros(size(V,2),size(V,3),size(V,4))];        
    else
        %Inertial acceleration
        inAcc = [x+2*V(2,:,:,:); y-2*V(1,:,:,:); zeros(size(V,2),size(V,3),size(V,4))];        
    end
    F = [V; inAcc-(mu1/R1.^3).*r1-(mu2/R2.^3).*r2];     %Time flow of the system
    
    %Update the differential configuration space vector
    dr = F;  
    
    %Reverse the flow for backward integration
    if (direction == -1)
        dr = -dr;
    end
end