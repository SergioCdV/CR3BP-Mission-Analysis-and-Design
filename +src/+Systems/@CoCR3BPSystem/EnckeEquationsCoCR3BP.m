%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: EnckeEquationsCoCR3BP.m 
% Issue: 0 
% Validated: 

%% Co-orbital CR3BP Dynamics in Encke form %%
% This function contains the description of the co-orbital CR3BP dynamics vector field

% Inputs: 
% Outputs: - vector ds, the differential vector field

% New versions: 

% Battin propagator for the three-body problem
function [ds] = EnckeEquationsCoCR3BP(t, j, s, u, params)
    % Define the initial phase space vector
    tgt = s(1:3,:);                         % Target synodic position vector
    r = s(7:9,:);                           % Relative synodic position vector
    x = r(1,:);                             % Relative synodic x coordinate
    y = r(2,:);                             % Relative synodic y coordinate 
    V = s(10:12,:);                         % Relative synodic velocity vector
    
    % Relevant system parameters
    mu = params(1);                         % Gravitational parameter of the system
    mup(1) = 1 - mu;                        % First primary normalized position
    mup(2) = mu;                            % Second primary normalized position
    R(:,1) = reshape(params(2:4), [], 1);   % Position vector of the first primary
    R(:,2) = reshape(params(5:7), [], 1);   % Position vector of the second primary

    % Inertial acceleration field
    gamma = [x; y; zeros(1,size(x,2))];                    
    gamma = gamma + [0 2 0; -2 0 0; 0 0 0] * V;
    ds = [V; gamma]; 

    % Gravitational forces, 
    for i = 1:length(mup)
        r_t = tgt - R(:,i);
        q = -dot((r + 2 * r_t), r, 1) ./ dot(r + r_t, r + r_t, 1);
        f = q .* (3 * (1 + q) + q.^2) ./ ( 1 + (1 + q).^(3/2) );

        % Encke acceleration method
        ds(4:6,:) = ds(4:6,:) - ( mup(i) / norm(r_t)^3 ) .* (f .* r_t + (1 + f) .* r);
    end

    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;
end