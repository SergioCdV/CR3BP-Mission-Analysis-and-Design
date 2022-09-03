%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/06/21
% File: rel_hamiltonian.m 
% Issue: 0 
% Validated: 

%% Relative Hamiltonian %%
% For a given gravitational parameter mu and a relative phase space vector, this function computes the relative
% Hamiltonain associated with that phase space vector

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - vector s, a 12x1 array containing both the position and velocity of the 
%           target mass and the relative particle respectively in the relative synodic frame

% Outputs: - scalar H, the relative Hamiltonian associated with the input phase space vector

% New versions: 

function [H] = rel_hamiltonian(mu, s)
    % Relative synodic position and velocity vectors
    x = s(7);                     % Relative synodic x coordinate
    y = s(8);                     % Relative synodic y coordinate
    v = s(10:12);                 % Relative velocity vector
    
    % Compute the kinetic energy 
    V = v-[y; x; 0];              % Total velocity vector
    T = (1/2)*dot(V,V);           % Kinetic energy
    
    % Compute the potential function
    U = rel_potential(mu, s);
    
    % Compute the relative Hamiltonian
    H = T+U;
end