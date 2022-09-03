%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/12/20
% File: inertial2synodic.m 
% Issue: 0 
% Validated: 

%% Inertial to Synodic %%
% For a given vector v at an epoch e, this function transforms a given vector 
% from the inertial frame to the synodic one or viceversa

% Inputs: - scalar epoch, epoch at which the vector r is known (in non-dimensional units)
%         - vector v, 3x1 array for position vector or 6x1 for a complete state vector
%         - boolean direction, 0 for inertial-to-synodic change and 1 for the viceversa conversion

% Outputs: - the vector V, v expressed in the new reference frame

% New versions: 

function [V] = inertial2synodic(epoch, v, direction)   
    % Constants 
    m = 6;                      % Phase space dimension
    omega = [0; 0; 1];          % Angular velocity of the synodic frame
    
    % Compute rotation matrix
    Q = [cos(epoch) sin(epoch) 0; -sin(epoch) cos(epoch) 0; 0 0 1];
    
    % Position transformation
    if (length(v) == m/2)
        if (direction)
            V = Q*v;            % Apply transformation
        else
            Q = Q.';                    
            V = Q*v;            % Apply transformation
        end
    
    %Velocity transformation
    elseif (length(v) == m)
        if (direction)
            V = v(4:6)-cross(omega, v(1:3));     % Synodic velocity vector 
            V = [Q*v(1:3); V];                   % Phase space vector
        else
            V = v(4:6)+cross(omega, v(1:3));     % Inertial velocity vector 
            V = [Q.'*v(1:3); V];                 % Phase space vector
        end
        
    %Error branch
    else
        error('No valid magnitude was selected');
    end    
end