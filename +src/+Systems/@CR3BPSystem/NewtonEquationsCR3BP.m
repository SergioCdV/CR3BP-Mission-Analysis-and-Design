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
    x = s(1,:);                       % Synodyc x coordinate
    y = s(2,:);                       % Synodyc y coordinate 
    z = s(3,:);                       % Synodyc z coordinate 
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
% 
%     if (0)
%         gamma = [x-2*V(2); y+2*V(1); 0];                                % Inertial acceleration
%     else
%         gamma = [x+2*V(2); y-2*V(1); 0];                                % Inertial acceleration
%     end
%     F = [V; gamma-(mup(1)/R(1)^3*r(:,1))-(mup(2)/R(2)^3*r(:,2))];       % Time flow of the system
    
    % Compute the GNC requirements 
%     if (~isempty(varargin))
%         if (~isempty(varargin{1}))
%             GNC = varargin{1};                      % GNC handling structure
%             if (iscell(GNC))
%                 GNC = GNC{1};
%             end
% 
%             % Include the GNC chain in the integration of the equations of motion
%             if (isfield(GNC.Algorithms, 'Control'))
%                 switch (GNC.Algorithms.Control)
%                     case 'MFKS'
%                         error('MFSK stationkeeping is not available for integration purposes')
%                     case 'HSK'
%                     otherwise
%                         error('No valid GNC algorithm was selected')
%                 end
%     
%                 % GNC scheme
%                 [~, ~, u] = GNCt_handler(GNC, s.', t);            % Compute the control law
%                 F(4:6) = F(4:6)+u;                                % Add the control vector 
%             end
%         end
%     end
    
    % Compute the variational equations if needed
%     if (flagVar)
%         % Compute the initial STM
%         Phi = reshape(s(n+1:end), [n n]);       % State Transition Matrix
%         J = abs_jacobian(mu,s);                 % Jacobian of the system 
%         dphi = J*Phi;                       	% Variational equations
%         dphi = reshape(dphi, [n^2 1]); 
%         
%         % Update the differential configuration space vector
%         dr = [F; dphi];
%     else
%         % Update the differential configuration space vector
%         dr = F;  
%     end
    
    % Reverse the flow for backward integration
%     if (direction == -1)
%         dr = -dr;
%     end
end