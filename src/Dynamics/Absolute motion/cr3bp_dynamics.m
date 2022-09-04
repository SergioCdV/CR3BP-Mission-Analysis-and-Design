%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: cr3bp_dynamics.m 
% Issue: 0 
% Validated: 

%% CR3BP Dynamics %%
% This function contains the 6-DOF dynamical model of the CR3BP dynamics 

% Inputs: - structure sysData, containing information about the central
%           bodies of the system, concerned mainly with dimensionless data.
%         - structure satData, containing information about the satellite
%           whose motion is beign studied: mass, moments of inertia, area dyadic...
%         - scalar t, an integration epoch. 
%         - vector flags, concerning selection of certain integration schemes.
%         - vector x, containing the n/n+n^2 degrees of freedom (phase space
%           and STM entries), depending on the type of integration.

% Output: - dx, an n/n+n^2 x 1 size vector containing the time flow of the full
%           dynamical model of the satellite.

% Methods: 

% New version: classical Euler solid rigid equations.

function [dx] = cr3bp_dynamics(sysData, satData, flags, t, x)
    %Define relevant constants 
    mu = sysData.mu;       %Reduced central bodies gravitational parameter
    attVar = flags(1);     %Boolean fo the integration of attitude dynamics
    stmVar = flags(2);     %Boolean for the integration of the traslational n+n^2 system
    
    %Integrate traslational motion
    dr = cr3bp_equations(mu, direction, stmVar, t, x);
    
    %Integrate attitude motion 
    if (attVar)
        dtheta = euler_equations(sysData, satData, t, x); 
    else
        dtheta = [];
    end
    
    %Update the differential state vector 
    dx = [dr; dtheta];
end