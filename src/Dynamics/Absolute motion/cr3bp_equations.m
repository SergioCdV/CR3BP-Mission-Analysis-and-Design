%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: cr3bp_equations.m 
% Issue: 0 
% Validated: 

%% CR3BP Dynamics %%
% This function contains the 42-DOF translational model used within the rest of
% the CR3BP General Library scripts. It accounts for a infinitesimal mass
% moving in the normalized, non dimensional synodic frame define by the two primaries, which
% are assumed to be in the same plane and in circular orbits. It also
% contains the integration of the first variational equations of the flow.

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - scalar direction (in binary format, 1 or -1), indicating the
%           time integration direction: 1 for forward integration, -1 for
%           backward integration.
%         - boolean flagVar, true for dyanmics and STM integration, 
%           false for only dynamical integration.
%         - scalar t, a reference epoch. 
%         - vector s, containing in an Nx1 array the phase space vector,
%           possibly augmented with the state transition matrix at time t.
%         - cell array varargin, to include GNC requirements on the motion
%           of the spacecraft

% Outputs: - vector dr, the differential vector field, which will include
%            the phase space trajectory and the STM integrated state when
%            flagVar is true.

% Methods: non-dimensional CR3BP dynamics in the synodic frame. 

% New versions: 

function [ds] = cr3bp_equations(mu, direction, flagVar, t, s, varargin)    
    %Equations of motion of the CR3BP
    global method_ID;

    switch (method_ID)
        %Deterministic models
        case 'Encke'
            ds = Encke_method(mu, direction, flagVar, t, s, varargin);    %Relative motion equations
        case 'Full nonlinear'
            ds = full_model(mu, direction, flagVar, t, s, varargin);      %Relative motion equations
        otherwise
            error('No valid model was chosen');
    end
end

%% Propagators
% Newtonian model
function [dr] = full_model(mu, direction, flagVar, t, s, varargin)
    %Constants 
    n = 6;                          %Phase space dimension 
    
    %Define the initial phase space vector
    x = s(1);                       %Synodyc x coordinate
    y = s(2);                       %Synodyc y coordinate 
    z = s(3);                       %Synodyc z coordinate 
    V = s(4:6);                     %Synodic velocity vector
    
    %Relevant system parameters
    mup(1) = 1-mu;                  %First primary normalized position
    mup(2) = mu;                    %Second primary normalized position
    r(:,1) = [(x+mup(2)); y; z];    %Relative position vector to the first primary
    r(:,2) = [(x-mup(1)); y; z];    %Relative position vector to the secondary primary
    R(1) = norm(r(:,1));            %Distance to the first primary
    R(2) = norm(r(:,2));            %Distance to the secondary primary
    
    %Compute the time flow of the system (depending on the time direction)
    if (direction == -1)
        gamma = [x-2*V(2); y+2*V(1); 0];                                %Inertial acceleration
    else
        gamma = [x+2*V(2); y-2*V(1); 0];                                %Inertial acceleration
    end
    F = [V; gamma-(mup(1)/R(1)^3*r(:,1))-(mup(2)/R(2)^3*r(:,2))];       %Time flow of the system
    
    %Compute the GNC requirements 
    GNC = varargin{1};                                                  %GNC handling structure
    if (~isempty(GNC))
        if (iscell(GNC))
            GNC = GNC{1};
        end

        %Include the GNC chain in the integration of the equations of motion
        if (isfield(GNC.Algorithms, 'Control'))
            switch (GNC.Algorithms.Control)
            end

            %GNC scheme
            [~, ~, u] = GNC_handler(GNC, s.', s.', t, true);            %Compute the control law
            F(4:6) = F(4:6)+u;                                          %Add the control vector 
        end
    end
    
    %Compute the variational equations if needed
    if (flagVar)
        %Compute the initial STM
        Phi = reshape(s(n+1:end), [n n]);       %State transition matrix
        J = abs_jacobian(mu,s);                 %Jacobian of the system 
        dphi = J*Phi;                       	%Variational equations
        dphi = reshape(dphi, [n^2 1]); 
        
        %Update the differential configuration space vector
        dr = [F; dphi];
    else
        %Update the differential configuration space vector
        dr = F;  
    end
    
    %Reverse the flow for backward integration
    if (direction == -1)
        dr = -dr;
    end
end

% Battin propagator for the three-body problem
function [dr] = Encke_method(mu, direction, flagVar, t, s, varargin)
    %Constants 
    n = 6;                          %Phase space dimension 
    L = libration_points(mu);       %Position array of the libration points
    L = L(1:3,2);                   %Position vector of the L1 point
    
    %Define the initial phase space vector
    r = s(1:3);                     %Synodic position vector
    V = s(4:6);                     %Synodic velocity vector
    x = r(1);                       %Synodyc x coordinate
    y = r(2);                       %Synodyc y coordinate 
    z = r(3);                       %Synodyc z coordinate 
    
    %Relevant system parameters
    mup(1) = 1-mu;                  %First primary normalized position
    mup(2) = mu;                    %Second primary normalized position
    R(:,1) = [-mu; 0; 0];           %Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];          %Synodic position of the second primary

    %Compute the time flow of the system (depending on the time direction)
    if (direction == -1)
        gamma = [x-2*V(2); y+2*V(1); 0];               %Inertial acceleration
    else
        gamma = [x+2*V(2); y-2*V(1); 0];               %Inertial acceleration
    end
    
    %Encke acceleration method
    for i = 1:length(mup)
        q = -dot((r+2*(L-R(:,i))),r)/norm(r+L-R(:,i))^2;
        f = q*(3*(1+q)+q^2)/(1+(1+q)^(3/2));
        gamma = gamma-mup(i)/norm(L-R(:,i))^3*(f*(L-R(:,i))+(1+f)*r);
    end
    F = [V; gamma];                         %Time flow of the system

    %Compute the GNC requirements 
    GNC = varargin{1};                      %GNC handling structure
    if (~isempty(GNC))
        if (iscell(GNC))
            GNC = GNC{1};
        end

        %Include the GNC chain in the integration of the equations of motion
        if (isfield(GNC.Algorithms, 'Control'))
            switch (GNC.Algorithms.Control)
            end

            %GNC scheme
            [~, ~, u] = GNC_handler(GNC, s.', s.', t, true);            %Compute the control law
            F(4:6) = F(4:6)+u;                                          %Add the control vector 
        end
    end
    
    %Compute the variational equations if needed
    if (flagVar)
        %Compute the initial STM
        Phi = reshape(s(n+1:end), [n n]);       %State transition matrix
        J = abs_jacobian(mu,s);                 %Jacobian of the system 
        dphi = J*Phi;                       	%Variational equations
        dphi = reshape(dphi, [n^2 1]); 
        
        %Update the differential configuration space vector
        dr = [F; dphi];
    else
        %Update the differential configuration space vector
        dr = F;  
    end
    
    %Reverse the flow for backward integration
    if (direction == -1)
        dr = -dr;
    end
end