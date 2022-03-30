%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/12/20
% File: synodic2lagrange.m 
% Issue: 0 
% Validated: 

%% Synodic to Lagrange %%
% For a given vector v at an epoch e, this function transforms a given vector 
% from the synodic frame to the Lagrange point reference frame.

% Inputs: - scalar mu, the reduced gravitational parameter of the system.
%         - vector v, a given position vector.
%         - scalar gamma, the characteristic distance of the libration point.
%         - boolean point, 1 for L1 or L3 and 2 for L2.
%         - boolean direction, 0 for synodic-to-lagrange change and 1 for the viceversa conversion.

% Outputs: - the vector V, v expressed in the new reference frame.

% New versions: 

function [V] = synodic2lagrange(mu, gamma, point, v, direction)
    %Libration point synodic position vector
    switch (point)
        case 1 
            X = [1-mu-gamma; 0; 0];
        case 2
            X = [1-mu+gamma; 0; 0];
        case 3
            X = [-mu+gamma; 0; 0];
        case 4 
            X = [mu+1/2; sqrt(3)/2; 0];
        case 5
            X = [mu+1/2; -sqrt(3)/2; 0];
    end
    
    %Main computations 
    if (direction == 0)
        V = (1/gamma)*(v-X);
    elseif (direction == 1)
        V = gamma*v+X;
    else
        error('No valid direction was selected');
    end
end