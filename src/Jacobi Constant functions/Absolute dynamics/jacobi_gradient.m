%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 21/02/21
% File: jacobi_gradient.m 
% Issue: 0 
% Validated: 

%% Jacobi gradient %%
% This function contains the algorithm to compute the vector gradient of the Jacobi constant for 
% a given system and state

% Inputs: - double mu, the reduced gravitational parameter of the system
%         - vector state, containing the system state variables at a given epoch

% Outputs: - vector T, the gradient of the Jacobi constant at a certain point

function [T] = jacobi_gradient(mu, s)
    % Constants of the system 
    mup(1) = 1-mu;                          % Reduced gravitational parameter of the first primary
    mup(2) = mu;                            % Reduced gravitational parameter of the second primary
    
    % State variables 
    x = s(1);                               % Synodic x coordinate
    y = s(2);                               % Synodic y coordinate
    z = s(3);                               % Synodic z coordinate
    v = s(4:end);                           % Synodic velocity vector
    
    % Main procedure 
    r(:,1) = [x+mup(2); y; z];                  % Relative position vector to the first primary
    r(:,2) = [x-mup(1); y; z];                  % Relative position vector to the second primary
    dU(:,1) = mup(1)/norm(r(:,1))^3*r(:,1);     % Acceleration due to the first primary
    dU(:,2) = mup(2)/norm(r(:,2))^3*r(:,2);     % Accleration due to the second primary
    dU = -[x; y; 0]+sum(dU,2);                  % Augmented potential gradient
    T = 2*[dU; v];                              % Jacobi Constant gradient vector 
end