%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: NewtonEquationsCoCR3BP.m 
% Issue: 0 
% Validated: 

%% Co-orbital CR3BP Dynamics %%
% This function contains the vector field of the co-orbital CR3BP system

% Inputs: 

% Outputs: - vector ds, the differential vector field of the system

% New versions: 

function [ds] = NewtonEquationsCoCR3BP(t, j, s, u, params)
    % Define the initial phase space vector
    r_t = s(1:3,:);                         % Target synodic position vector
    r_r = s(7:9,:);                         % Relative synodic position vector
    x = r_r(1,:);                           % Relative synodic x coordinate
    y = r_r(2,:);                           % Relative synodic y coordinate 
    V = s(10:12,:);                         % Relative synodic velocity vector
    
    % Relevant system parameters
    mu = params(1);                                    % Gravitational parameter of the system
    mup(1) = 1 - mu;                                   % First primary normalized position
    mup(2) = mu;                                       % Second primary normalized position
    R(:,1) = reshape(params(2:4), [], 1);              % Position vector of the first primary
    R(:,2) = reshape(params(5:7), [], 1);              % Position vector of the second primary

    Rr(1:3,:) = r_t - R(:,1);                          % Synodic relative position of the target to the first primary
    Rr(4:6,:) = r_t - R(:,2);                          % Synodic relative position of the target to the second primary

    Rr_norm(1,:) = sqrt( dot(Rr(1:3,:), Rr(1:3,:), 1) );     % Distance to the first primary
    Rr_norm(2,:) = sqrt( dot(Rr(4:6,:), Rr(4:6,:), 1) );     % Distance to the secondary primary
    
    % Compute the time flow of the system
    gamma = [x; y; zeros(1,size(x,2))];                % Inertial acceleration terms
    gamma = gamma + [0 2 0; -2 0 0; 0 0 0] * V;
    ds = [V; gamma]; 
 
    % Gravitational forces
    F =   + mup(1) * ( Rr(1:3,:) ./ Rr_norm(1,:).^3 - (r_r + Rr(1:3,:)) ./ sqrt( dot(r_r + Rr(1:3,:), r_r + Rr(1:3,:), 1) ).^3 );
    F = F + mup(2) * ( Rr(4:6,:) ./ Rr_norm(2,:).^3 - (r_r + Rr(4:6,:)) ./ sqrt( dot(r_r + Rr(4:6,:), r_r + Rr(4:6,:), 1) ).^3 );
    ds(4:6,:) = ds(4:6,:) + F;
 
    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;
end