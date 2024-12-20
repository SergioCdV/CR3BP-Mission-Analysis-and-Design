%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 19/12/24
% File: AugmentedPotential.m 
% Issue: 0 
% Validated: 

%% Augmented potential %%
% For a given gravitational parameter mu and position vectors r, this function computes the 
% augmented potential function associated with that input position vector

% Inputs: - double mu, the reduced gravitational parameter of the system 
%         - array r [3xN], containing the synodic position vectors

% Outputs: - vector U [N], the augmented potential function associated with the input position vector 

% New versions:

function [U] = AugmentedPotential(mu, r)
    % Constants of the problem 
    mup(1) = 1 - mu;                    % Gravitational parameters of the first primary
    mup(2) = mu;                        % Gravitational parameters of the second primary
    
    % Obtain synodic coordinates 
    x = r(1,:);                         % Synodic x coordinate
    y = r(2,:);                         % Synodic y coordinate 
    z = r(3,:);                         % Synodic z coordinate
    
    % Compute distance to the primaries
    R(1:3,:) = [x + mup(2); y; z];      % Relative position vector to the first primary    
    R(4:6,:) = [x - mup(1); y; z];      % Relative position vector to the second primary 
    
    % Augmented potential function 
    U = -0.5 * (x.^2 + y.^2) - ( mup(1) ./ sqrt(dot(R(1:3,:), R(1:3,:), 1)) ) - ( mup(2) ./ sqrt(dot(R(4:6,:), R(4:6,:), 1)) );

    U = U - 0.5 * mup(1) * mup(2);
end