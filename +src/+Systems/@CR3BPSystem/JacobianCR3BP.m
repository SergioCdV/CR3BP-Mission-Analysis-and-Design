%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 23/12/24
% File: JacobianCR3BP.m 
% Issue: 0 
% Validated: 

%% Jacobian for the CR3BP %%
% This function contains the jacobian matrix of the vector field of the CR3BP. It accounts for a infinitesimal mass
% moving in the normalized, non dimensional synodic frame define by the two primaries, which
% are assumed to be in the same plane and in circular orbits

% Inputs: - scalar mu, the reduced gravitational parameter of the system 
%         - vector s, containing in an Mx1 array the phase space vector, with M 
%           the phase space dimension

% Outputs: - matrix J, containing the Jacobian evaluated at the phase space vector

% Methods: non-dimensional CR3BP dynamics in the synodic frame 

% New versions: 

function [J] = JacobianCR3BP(mu, s) 
    % Define the phase space vector
    x = s(1,:);                       % Synodyc x coordinate
    y = s(2,:);                       % Synodyc y coordinate 
    z = s(3,:);                       % Synodyc z coordinate 
    
    % Relevant system parameters
    mup(1) = 1 - mu;                  % First primary normalized position
    mup(2) = mu;                      % Second primary normalized position

    r(1:3,:) = [x + mup(2); y; z];    % Relative position vector to the first primary
    r(4:6,:) = [x - mup(1); y; z];    % Relative position vector to the secondary primary

    R(1,:) = sqrt( dot(r(1:3,:), r(1:3,:), 1) );            % Distance to the first primary
    R(2,:) = sqrt( dot(r(4:6,:), r(4:6,:), 1) );            % Distance to the secondary primary
    
    % Constants
    O = zeros(3,3);                    % Null matrix
    I = eye(3);                        % Identity matrix
    K = [0 2 0; -2 0 0; 0 0 0];        % Coriolis dyadic

    % Pre-allocation 
    J = zeros(6, 6 * size(s,2));

    % First variations of the augmented potential function (Hessian of the potential)
    G(1,:) = 1 - ( mup(1) ./ R(1,:).^3) .* ( 1 - 3 .* ( r(1,:) ./ R(1,:) ).^2 ) - ( mup(2) ./ R(2,:).^3 ) .* (1 - 3 .* ( r(4,:) ./ R(2,:) ).^2 ); 
    G(2,:) = 3 * r(2,:) .* ( (mup(1) ./ R(1,:).^5) .* r(1,:) + (mup(2) ./ R(2,:).^5) .* r(4,:) );
    G(3,:) = 3 * r(3,:) .* ( (mup(1) ./ R(1,:).^5) .* r(1,:) + (mup(2) ./ R(2,:).^5) .* r(4,:) );

    G(4,:) = G(2,:);
    G(5,:) = 1 - ( mup(1) ./ R(1,:).^3) .* ( 1 - 3 .* ( r(2,:) ./ R(1,:) ).^2 ) - ( mup(2) ./ R(2,:).^3 ) .* (1 - 3 .* ( r(5,:) ./ R(2,:) ).^2 );
    G(6,:) = 3 * r(2,:) .* ( ( mup(1) ./ R(1,:).^5 ) .* r(3,:) + ( mup(2) ./ R(2,:).^5 ) .* r(3,:) );

    G(7,:) = G(3,:); 
    G(8,:) = G(6,:);
    G(9,:) = 0 - ( mup(1) ./ R(1,:).^3) .* ( 1 - 3 .* ( r(3,:) ./ R(1,:) ).^2 ) - ( mup(2) ./ R(2,:).^3 ) .* (1 - 3 .* ( r(6,:) ./ R(2,:) ).^2 );

    % Compute the first variational equations evaluated at the reference
    for i = 1:size(s,2)
        % Jacobian of the system 
        idx = 1 + 6 * (i - 1) : 6 * i;
        J(:,idx) = [O I; reshape(G(:,i), 3, 3) K];             
    end
end