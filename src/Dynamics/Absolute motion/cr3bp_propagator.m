%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/10/20
% File: cr3bp_propagator.m 
% Issue: 0 
% Validated: 

%% CR3BP Propagator %%
% This function contains the propagator to integrate the CR3BP dynamics using different methods and numerical schemes.

% Inputs: - structure setup, selecting the type of integration to be
%           performed.
%         - scalar mu, the reduced gravitational parameter of the system. 
%         - scalar direction (in binary format, 1 or -1), indicating the
%           time integration direction: 1 for forward integration, -1 for
%           backward integration.
%         - boolean flagVar, true for dyanmics and STM integration, 
%           false for only dynamical integration.
%         - scalar t, a reference epoch. 
%         - vector s, containing in an Nx1 array the phase space vector,
%           possibly augmented with the state transition matrix at time t.
%         - cell array varargin, to include GNC requirements on the motion
%           of the spacecraft

% Outputs: - vector dr, the differential vector field, which will include
%            the phase space trajectory and the STM integrated state when
%            flagVar is true.

% Methods: non-dimensional CR3BP dynamics in the synodic frame. 

% New versions: 

function [ds] = cr3bp_propagator(setup, mu, direction, flagVar, t, s, varargin)    
    %Equations of motion of the CR3BP
    method_ID = setup.Method;

    switch (method_ID)
        %Deterministic models
        case 'Newton'    
            ds = cr3bp_equations(mu, direction, flagVar, t, s, varargin);            %Absolute equations of motion
            
        case 'Encke'
            L = libration_points(mu);                                                %System libration points (varargin)
            d = repmat(s(1:3),1,size(L,2))-L(1:3,:);                                 %Relative distance to the points
            D = sqrt(dot(d,d,1));                                                    %Relative distance to the libration point
            L = L(1:3, D == min(D));                                                 %Lagrange point to which relative motion is computed (minimum distance)
            s(1:3) = d(:,D == min(D));                                               %Relative distance to the libration point
            ds = encke_dynamics(mu, L, direction, flagVar, t, s, varargin);          %Relative motion equations

        otherwise
            error('No valid integration scheme was chosen.');
    end
end