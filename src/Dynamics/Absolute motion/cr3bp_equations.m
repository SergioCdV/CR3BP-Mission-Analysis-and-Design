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

function [dr] = cr3bp_equations(mu, direction, flagVar, t, s, varargin)
    %Constants 
    n = 6;       %Phase space dimension 
    
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
    if (~isempty(varargin))
        GNC = varargin{1};                                              %GNC handling structure
        if (iscell(GNC))
            GNC = GNC{1};
        end

        %Integrate the relative position for the PID controller
        if (isfield(GNC.Algorithms, 'Control'))
            switch (GNC.Algorithms.Control)
            end

            %GNC scheme
            [~, ~, u] = GNC_handler(GNC, s.', s.', t, true);            %Compute the control law
            F(4:6) = F(4:6)+u;                                          %Add the control vector 
        end

        %Change to a more numerically stable propagator
        if (isfield(GNC.Algorithms, 'Propagation'))
            switch (GNC.Algorithms.Propagation)
                case 'Battin'
                    F = Battin_method(mu, s);
                otherwise
            end
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

%% Auxiliary propagators 
% Battin propagator for the three-body problem
function [drho] = Battin_method(mu, s_t)
    %Constants of the system 
    mu_r(1) = 1-mu;               %Reduced gravitational parameter of the first primary 
    mu_r(2) = mu;                 %Reduced gravitational parameter of the second primary 
    
    %State variables 
    r_t = s_t(1:3);               %Synodic relative position 
    v_t = s_t(4:6);               %Synodic relative velocity 
    
    %Synodic position of the primaries 
    R(1:3,1) = [-mu; 0; 0];       %Synodic position of the first primary
    R(1:3,2) = [1-mu; 0; 0];      %Synodic position of the second primary
    
    %Encke acceleration method
    gamma = [2*v_t(2)+r_t(1); -2*v_t(1)+r_t(2); 0];
    for i = 1:length(mu_r)
        q = -dot(2*(-R(:,i))+r_t,r_t)/norm(r_t-R(:,i))^2;
        f = q*(3*(1+q)+q^2)/(1+(1+q)^(3/2));
        gamma = gamma - (mu_r(i)/norm(R(:,i))^3)*(f*(-R(:,i))+(1+f)*r_t);
    end
    
    %Equations of motion 
    drho = [v_t; 
            gamma];
end