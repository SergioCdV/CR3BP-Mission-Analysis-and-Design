%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 19/12/24
% File: JacobiGradient.m 
% Issue: 0 
% Validated: 

%% Jacobi Gradient %%
% This function contains the algorithm to compute the vector gradient of the Jacobi constant for 
% a given system and state

% Inputs: - double mu, the reduced gravitational parameter of the system
%         - vector state [6xN], containing the system state variables at a given epoch

% Outputs: - vector T[6xN], the gradient of the Jacobi constant at a certain point

function [T] = JacobiGradient(mu, s)
    % Constants of the system 
    mup(1) = 1 - mu;                            % Reduced gravitational parameter of the first primary
    mup(2) = mu;                                % Reduced gravitational parameter of the second primary
    
    % State variables 
    x = s(1,:);                                 % Synodic x coordinate
    y = s(2,:);                                 % Synodic y coordinate
    z = s(3,:);                                 % Synodic z coordinate
    v = s(4:6,:);                               % Synodic velocity vector [vx, vy, vz]
    
    % Main procedure 
    r(1:3,:) = [x + mup(2); y; z];              % Relative position vector to the first primary
    r(4:6,:) = [x - mup(1); y; z];              % Relative position vector to the second primary

    dU(1:3,:) = mup(1) ./ sqrt(dot(r(4:6,:), r(4:6,:), 1)).^3 .* r(4:6,:);     % Acceleration due to the first primary
    dU(4:6,:) = mup(2) ./ sqrt(dot(r(4:6,:), r(4:6,:), 1)).^3 .* r(4:6,:);     % Accleration due to the second primary
    
    dU = -[x; y; 0] + dU(1:3,:) + dU(4:6,:);    % Augmented potential gradient
    T = 2 * [dU; v];                            % Jacobi Constant gradient vector 
end