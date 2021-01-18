%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/09/20
% File: jacobi_constant.m 
% Issue: 0 
% Validated: 

%% Jacobi constant %%
% For a given gravitational parameter mu, a particular vector velocity v
% and a position vector r, this function computes the Jacobi constant
% associated with that phase space vector. 

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - vector s, a 6x1 array containing both position and velocity of the satellite 
%           respectively in the synodic frame.

% Outputs: - scalar J, the Jacobi Constant associated with the input phase space vector. 

% New versions: 

function [J] = jacobi_constant(mu, s)
    %Define the synodic position and velocity vectors
    r = s(1:3);     %Position vector
    v = s(4:6);     %Velocity vector
        
    %Compute the augmented potential function
    U = augmented_potential(mu, r);
    
    %Compute the Jacobi Constant
    J = 2*U-dot(v,v);
end