%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 23/12/24
% File: NewtonEquationsCR3BP.m 
% Issue: 0 
% Validated: 

%% CR3BP Dynamics %%
% This function contains the vector field of the CR3BP system. It accounts for a infinitesimal mass
% moving in the normalized, non dimensional synodic frame define by the two primaries, which
% are assumed to be in the same plane and in circular orbits. It also
% contains the integration of the first variational equations of the flow

% Inputs: 

% Outputs: - vector ds, the differential vector field of the system

% New versions: 

function [ds] = NewtonEquationsCR3BP(t, j, s, u, params)
    % Define the initial phase space vector
    x = s(1,:);                       % Synodic x coordinate
    y = s(2,:);                       % Synodic y coordinate 
    z = s(3,:);                       % Synodic z coordinate 
    V = s(4:6,:);                     % Synodic velocity vector
    
    % Relevant system parameters
    mu = params(1);                   % Gravitational parameter of the system
    mup(1) = 1 - mu;                  % First primary normalized position
    mup(2) = mu;                      % Second primary normalized position

    r(1:3,:) = [x + mup(2); y; z];    % Relative position vector to the first primary
    r(4:6,:) = [x - mup(1); y; z];    % Relative position vector to the secondary primary

    R(1,:) = sqrt( dot(r(1:3,:), r(1:3,:), 1) );            % Distance to the first primary
    R(2,:) = sqrt( dot(r(4:6,:), r(4:6,:), 1) );            % Distance to the secondary primary
    
    % Compute the time flow of the system
    gamma = [x; y; zeros(1,size(x,2))];                     % Inertial acceleration terms
    gamma = gamma + [0 2 0; -2 0 0; 0 0 0] * V;
    ds = [V; gamma]; 

    % Gravitational forces
    ds(4:6,:) = ds(4:6,:) - mup(1) ./ R(1,:).^3 .* r(1:3,:) - mup(2) ./ R(2,:).^3 .* r(4:6,:);

    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;
end