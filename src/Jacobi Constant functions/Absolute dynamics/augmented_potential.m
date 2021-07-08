%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/09/20
% File: augmented_potential.m 
% Issue: 0 
% Validated: 

%% Jacobi constant %%
% For a given gravitational parameter mu and position vector r, this function computes the 
% potential function associated with that input position vector. 

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - vector r, a 3x1 array containing the synodic position vector.

% Outputs: - scalar U, the augmented potential function associated with the input position vector. 

% New versions:

function [U] = augmented_potential(mu, r)
    %Constants of the problem 
    mup(1) = 1-mu;                      %Gravitational parameters of the most massive primary
    mup(2) = mu;                        %Gravitational parameters of the least massive primary 
    
    %Obtain synodic coordinates 
    x = r(1);                           %Synodic x coordinate
    y = r(2);                           %Synodic y coordinate 
    z = r(3);                           %Synodic z coordinate
    
    %Compute distance to the primaries
    r(:,1) = [x+mup(2); y; z];          %Position vector to the most massive primary    
    r(:,2) = [x-mup(1); y; z];          %Position vector to the least massive primary
    
    %Augmented potential function 
    U = -(1/2)*(x^2+y^2)-(mup(1)/norm(r(:,1)))-(mup(2)/norm(r(:,2)));
end