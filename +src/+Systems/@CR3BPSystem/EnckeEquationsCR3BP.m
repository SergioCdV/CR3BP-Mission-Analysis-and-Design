%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 23/12/24
% File: EnckeEquationsCR3BP.m 
% Issue: 0 
% Validated: 

%% CR3BP Dynamics in Encke form %%
% This function contains the description of the CR3BP dynamics vector field. It accounts for a infinitesimal mass
% moving in the normalized, non dimensional synodic frame define by the two primaries, which
% are assumed to be in the same plane and in circular orbits. It also
% contains the integration of the first variational equations of the flow

% Inputs: 
% Outputs: - vector ds, the differential vector field

% New versions: 

% Battin propagator for the three-body problem
function [ds] = EnckeEquationsCR3BP(t, j, s, u, params)
    % Define the initial phase space vector
    L = s(1:3,:);                           % Relative position of the target 
    s = s(4:end,:);                         % Relative state vector
    r = s(1:3,:);                           % Synodic position vector
    V = s(4:6,:);                           % Synodic velocity vector
    x = r(1,:);                             % Synodic x coordinate
    y = r(2,:);                             % Synodic y coordinate 
    
    % Relevant system parameters
    mu = params(1);                         % Gravitational parameter of the system
    mup(1) = 1 - mu;                        % First primary normalized position
    mup(2) = mu;                            % Second primary normalized position
    R(:,1) = reshape(params(2:4), [], 1);   % Position vector of the first primary
    R(:,2) = reshape(params(2:4), [], 1);   % Position vector of the second primary

    % Inertial acceleration field
    gamma = [x; y; zeros(1,size(x,2))];                    
    gamma = gamma + [0 2 0; -2 0 0; 0 0 0] * V;
    ds = [V; gamma]; 

    % Gravitational forces, 
    for i = 1:length(mup)
        r_t = L - R(:,i);
        q = -dot((r + 2 * r_t), r, 1) ./ dot( r + r_t, r + r_t, 1);
        f = q .* (3 * (1 + q) + q.^2) ./ ( 1 + (1 + q).^(3/2) );

        % Encke acceleration method
        ds(4:6,:) = ds(4:6,:) - ( mup(i) / norm(r_t)^3 ) .* (f .* r_t + (1 + f) .* r);
    end

    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;

%     % Compute the GNC requirements    
%     if (~isempty(varargin))
%         if (~isempty(varargin{1}))
%             GNC = varargin{1};              % GNC handling structure
%             if (iscell(GNC))
%                 GNC = GNC{1};
%             end
% 
%             % Include the GNC chain in the integration of the equations of motion
%             if (isfield(GNC.Algorithms, 'Control'))
%                 switch (GNC.Algorithms.Control)
%                 end
%     
%                 % GNC scheme
%                 [~, ~, u] = GNC_handler(GNC, s.', s.', t, true);            % Compute the control law
%                 F(4:6) = F(4:6)+u;                                          % Add the control vector 
%             end
%         end
%     end
%     
%     % Compute the variational equations if needed
%     if (flagVar)
%         % Compute the initial STM
%         Phi = reshape(s(n+1:end), [n n]);       % State transition matrix
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
%     
%     % Reverse the flow for backward integration
%     if (direction == -1)
%         dr = -dr;
%     end
end