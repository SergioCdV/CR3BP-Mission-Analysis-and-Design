%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 21/02/21
% File: jacobi_gradient.m 
% Issue: 0 
% Validated: 

%% Jacobi gradient %%
% This function contains the algorithm to compute the vector gradient of the Jacobi constant for 
% a given system and state.

% Inputs: - double mu, the reduced gravitational parameter of the system.
%         - vector state, containing the system state variables at a given epoch. 

% Outputs: - vector T, the gradient of the Jacobi constant at a certain point.

function [T] = jacobi_gradient(mu, state)
    %Constants of the system 
    mu1 = 1-mu;     %Reduced gravitational parameter of the first primary
    mu2 = mu;       %Reduced gravitational parameter of the secnd primary
    
    %State variables 
    x = state(1);       %Synodic x coordinate
    y = state(2);       %Synodic y coordinate
    z = state(3);       %Synodic z coordinate
    v = state(4:end);   %Synodic velocity vector
    
    %Main procedure 
    r1 = [(x+mu2); y; z];               %Relative position vector to the first primary
    r2 = [(x-mu1); y; z];               %Relative position vector to the secondary primary
    dU1 = -mu1/norm(r1)^3*r1;           %Acceleration due to the first primary
    dU2 = -mu2/norm(r2)^3*r2;           %Accleration due to the second primary
    dU = [x; y; 0]-dU1-dU2;             %Augmented potential gradient
    T = 2*[dU; v];                       %Jacobi Constant gradient vector 
end