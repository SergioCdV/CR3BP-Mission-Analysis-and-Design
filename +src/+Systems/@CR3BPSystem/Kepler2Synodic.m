%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 22/12/24
% File: Kepler2Synodic.m 
% Issue: 0 
% Validated: 

%% Synodic to Inertial %%
% For a given CR3BP, this function computes the homogeneous transformation
% (4x4) to transform from the synodic to the inertial reference frame

% Inputs: - scalar mu, the gravitational parameter of the system
%         - scalar idx, the ID of the primary of interest
%         - vector theta, 1xm, the variable describing the motion of the second
%           primary with respect to the primary
%         - boolean direction, denoting the direction of the transformation
%           (true for synodic to inertial)

% Outputs: - matrix T, the 4x4 x m homogeneous matrix defining the
%            transformation (Galilean group)

% New versions: 

function [T] = Kepler2Synodic(mu, idx, theta, direction)
    % Preallocation 
    T = src.Systems.CR3BPSystem.Synodic2Inertial(theta, direction);

    if (idx == 1)
        R = -mu; 
    else
        R = 1 - mu;
    end

    if ( direction )
        for i = 1:length(theta)
            T(1:3, 4*i) = [+R; 0; 0];
        end
    else
        for i = 1:length(theta)
            T(1:3, 4*i) = [-R; 0; 0];
        end
    end
end
