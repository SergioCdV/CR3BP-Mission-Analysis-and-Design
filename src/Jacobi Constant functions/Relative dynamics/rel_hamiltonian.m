%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/06/21
% File: rel_hamiltonian.m 
% Issue: 0 
% Validated: 

%% Relative Hamiltonian %%
% For a given gravitational parameter mu, a particular relative velocity
% and a position vectors v and r, this function computes the relative
% Hamiltonain associated with that phase space vector. 

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - vector s, a 12x1 array containing both the position and velocity of the 
%           target mass and the relative particle respectively in the relative synodic frame.

% Outputs: - scalar H, the relative Hamiltonian associated with the input phase space vector. 

% New versions: 

function [H] = rel_hamiltonian(mu, s)
    %Define the relative synodic position and velocity vectors
    v = s(12:15);                 %Relative velocity vector
    
    %Compute the kinetic energy 
    V = v-[x; y; 0];              %Total velocity vector
    T = (1/2)*norm(V)^2;          %Kinetic energy
    
    %Compute the augmented potential function
    U = rel_potential(mu, s);
    
    %Compute the relative Hamiltonian
    H = T+U;
end