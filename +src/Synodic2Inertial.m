%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 17/06/24
% File: SynodicTransformation.m 
% Issue: 0 
% Validated: 

%% Synodic to Inertial %%
% For a given CR3BP, this function computes the homogeneous transformation
% (4x4) to transform from the synodic to the inertial reference frame

% Inputs: - scalar mu, the parameter of the system
%         - vector theta, 1xm, the variable describing the motion of the second
%           primary with respect to the primary
%         - boolean direction, denoting the direction of the transformation
%           (true for synodic to inertial)

% Outputs: - matrix T, the 4x4 x m homogeneous matrix defining the
%            transformation (Galilean group)

% New versions: 

function [T] = Synodic2Inertial(mu, theta, direction)
    % Preallocation 
    T = zeros(4, 4 * size(theta,2));

    % Compute the transformation from the synodic to the inertial reference frame
    for i = 1:size(theta,2)
        % Compute the rotation matrix 
        cos_theta = cos( theta(i) ); 
        sin_theta = sin( theta(i) );

        T(1:3,1+3*(i-1):3*i) = [cos_theta -sin_theta 0; sin_theta cos_theta 0; 0 0 1];

        % Compute the displacement
        T(4,4 * i) = 1;

        % Consider the direction of the transformation
        if ~direction
            T(1:3, 1+3*(i-1):3*i) = T(1:3,1+3*(i-1):3*i).';      % Inverse of the rotation matrix
        end
    end
end
