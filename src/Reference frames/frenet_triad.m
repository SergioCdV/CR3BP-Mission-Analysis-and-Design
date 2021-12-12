%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 13/03/21
% File: frenet_triad.m 
% Issue: 0 
% Validated: 

%% Frenet triad %%
% For a given gravitational parameter mu and a space state vector, for the CR3BP dynamics, this function
% computes the Frenet-Serret triad evaluated at that phase point.

% Inputs: - scalar mu, the reduced gravitational parameter of the system.
%         - vector s, phase space vector of at minimum 6x1.

% Outputs: - the array T, containing the three unitary vectors of the Frenet triad. It is also
%            the rotation matrix from the synodic frame to the F-S frame.

% New versions: .

function [T] = frenet_triad(mu, s)
    %Evaluate the vector field at the phase point
    [f] = cr3bp_equations(mu, 1, false, 0, s.');    %Vector field
    
    %Compute the Frenet-Serret triad 
    t = f(1:3)/norm(f(1:3));                        %Tangent vector 
    n = cross(f(1:3),cross(f(4:6), f(1:3)));        %Normal vecor
    n = n/norm(n);                                  %Normal vector
    b = cross(t, n);                                %Binormal vector
    
    %Output 
    T = [t n b];
end