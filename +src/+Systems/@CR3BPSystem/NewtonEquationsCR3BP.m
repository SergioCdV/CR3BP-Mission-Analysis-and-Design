%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 23/12/24
% File: NewtonEquationsCR3BP.m 
% Issue: 0 
% Validated: 

%% CR3BP Dynamics %%

% Inputs: 
% Outputs: - vector ds, the differential vector field of the system

% New versions: 

function [ds] = NewtonEquationsCR3BP(t, j, s, u, params)
    % Define the initial phase space vector
    r_t = s(1:3,:);                   % Synodic position vector
    x = s(1,:);                       % Synodic x coordinate
    y = s(2,:);                       % Synodic y coordinate 
    z = s(3,:);                       % Synodic z coordinate 
    V = s(4:6,:);                     % Synodic velocity vector
    
    % Relevant system parameters
    mu = params(1);                                    % Gravitational parameter of the system
    mup(1) = 1 - mu;                                   % First primary normalized position
    mup(2) = mu;                                       % Second primary normalized position
    Rp(:,1) = reshape(params(2:4), [], 1);             % Position vector of the first primary
    Rp(:,2) = reshape(params(5:7), [], 1);             % Position vector of the second primary
    
    r(1:3,:) = r_t(1:3,:) - Rp(:,1);                   % Synodic relative position of the target to the first primary
    r(4:6,:) = r_t(1:3,:) - Rp(:,2);                   % Synodic relative position of the target to the second primary

    R(1,:) = sqrt( dot(r(1:3,:), r(1:3,:), 1) );       % Distance to the first primary
    R(2,:) = sqrt( dot(r(4:6,:), r(4:6,:), 1) );       % Distance to the secondary primary
    
    % Compute the time flow of the system
    gamma = [x; y; zeros(1,size(x,2))];                % Inertial acceleration terms
    gamma = gamma + [0 2 0; -2 0 0; 0 0 0] * V;
    ds = [V; gamma]; 

    % Gravitational forces
    Accg = - mup(1) ./ R(1,:).^3 .* r(1:3,:) - mup(2) ./ R(2,:).^3 .* r(4:6,:);
    ds(4:6,:) = ds(4:6,:) + Accg;

    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;
end