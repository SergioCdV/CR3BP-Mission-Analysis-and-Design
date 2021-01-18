%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/10/20
% File: timeflow.m 
% Issue: 0 
% Validated: 

%% Time flow %%
% This function contains a routine to obtain the flow of a initial condition from t0 to time 
% t1. It also outputs the STM at that point (stroboscopic map).

% Input: - scalar t0, an initial time from which integration is applied. 
%        - scalar tf, a final time to end integration.
%        - scalar mu, reduced gravitational parameter of the system.
%        - scalar direction, corcening the time direction of the
%          integration: 1 for forward integration, -1 for backward
%          integration.
%        - vector x, with sizes 42x1, containing initial conditions for the
%          integration.

% Outputs: - vector xf, phase space vector at time tf  (6x1).
%          - matrix STM, the state transition matrix of the system at time
%          tf.

% New versions: integration with attitude dynamics

function [xf, STM] = timeflow(t0, tf, direction, mu, x) 
    %Set general integration parameters 
    RelTol = 2.5e-14; 
    AbsTol = 1e-22; 
    options = odeset('RelTol', RelTol, 'AbsTol', AbsTol); 
    
    %Constants 
    n = 6;                                  %Phase space dimension 
        
    %Time integration vector 
    rho = 5000;                             %Number of time steps
    if (direction ~= -1)
        tspan = linspace(0, (tf-t0), rho);
    else
        tspan = linspace(0, (t0-tf), rho);
    end
    
    %Dynamics integration 
    [~, xf] = ode113(@(t,s)cr3bp_equations(mu, direction, true, t, s), tspan, x, options);
    
    %State transition matrix         
    STM = reshape(shiftdim(xf(end,n+1:end)), n, n);
end