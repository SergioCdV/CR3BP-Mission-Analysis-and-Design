%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 19/12/24
% File: JacobiConstant.m 
% Issue: 0 
% Validated: 

%% Jacobi constant %%
% For a given gravitational parameter mu, phase space vector s, this function computes the Jacobi constant
% associated with that phase space vector

% Inputs: - double mu, the reduced gravitational parameter of the system.
%         - array s [6xN] containing N states, composed of both positions and velocities in the synodic frame

% Outputs: - vector J [N], the Jacobi Constant associated with the input phase space vector
%          - vector H [N], the Hamiltonian associated with the input phase space vector

function [J, H] = JacobiConstant(mu, s)
    % Define the synodic position and velocity vectors
    r = s(1:3,:);     % Position vector
    v = s(4:6,:);     % Velocity vector
        
    % Compute the augmented potential function
    U = src.Systems.CR3BPSystem.AugmentedPotential(mu, r);
    
    % Compute the Jacobi Constant of the input vectors
    J = -2 * U - dot(v, v, 1);
    H = - 0.5 * J;
end