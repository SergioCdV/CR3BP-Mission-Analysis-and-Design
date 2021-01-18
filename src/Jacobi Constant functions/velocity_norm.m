%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/09/20
% File: velocity_norm.m 
% Issue: 0 
% Validated: 

%% Velocity norm %%
% For a given gravitational parameter mu, the velocity_norm function
% computes the magnitude of the velocity vector associated with a particular Jacobi Constant 
% and a specified position vector.

% Inputs: - scalar mu, the reduced gravitational parameter of the system.
%         - vector r, a position vector expressed in the non-dimensional
%           synodic frame.
%         - scalar J, the associated Jacobi constant.

% Outputs: - the scalar V, the magnitude of the velocity vector.

% New versions: 

function [V] = velocity_norm(mu, r, J) 
    %Compute the augmented potential function
    U = augmented_potential(mu, r); 
    
    %Compute the magnitude of the synodic velocity
    V = sqrt(2*U-J);
end