%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 10/03/21
% File: fulrel_motion.m 
% Issue: 0 
% Validated: 

%% Full relative motion in the CR3BP %%
% This function contains the full nonlinear relative motion between two punctual particles in the 
% CR3BP. It accounts for the two masses moving in the normalized, non dimensional synodic frame 
% define by the two primaries, which are assumed to be in the same plane and in circular orbits

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - scalar direction (in binary format, 1 or -1), indicating the
%           time integration direction: 1 for forward integration, -1 for
%           backward integration
%         - boolean relFlagVar, true for relative dyanmics and STM integration, 
%           false for only dynamical integration
%         - scalar omega, the angular velocity of the target torus
%         - function handle g, the torus function of the target periodic
%           orbit
%         - scalar theta, a reference torus angle
%         - vector s, containing in an Nx1 array the phase space vector,
%           possibly augmented with the state transition matrix at time t
%         - cell array varargin, to include GNC requirements on the motion
%           of the spacecraft

% Outputs: - vector ds, the differential vector field, which will include
%            the phase space trajectory

% Methods: 

% New versions: 

function [drho] = geo_model(mu, direction, relFlagVar, omega, g, theta, s, varargin)
    % Equations of motion of the target
    s_t = feval(g, theta);                          % Target state
    
    % Equations of motion of the relative state 
    drho = geometrical_model(mu, omega, s_t, s);    % Relative motion equations

    % Augmented integral state dynamics
    if (mod(length(s), 9) == 0)
        m = 9;                                      % Augmented phase space dimension
        drho = [drho; s_r(7:9)];                    % Integrate the relative position for the PID controller
    end
    
    % GNC handler 
    if (~isempty(varargin))
        GNC = varargin{1};                                        % GNC handling structure
        if (iscell(GNC))
            GNC = GNC{1};
        end
    else
        GNC = [];                                                 % Empty GNC structure
    end

    if (~isempty(GNC))
        % GNC scheme
        [~, ~, u] = GNC_handler(GNC, s_t(1:6).', s.', t);         % Compute the control law
        drho(4:6) = drho(4:6)+u;                                  % Add the control vector       
    end
    
    % Relative variational equations
    if (relFlagVar)
        dphi = var_equations(mu, t, s);      % Linear variational equations
        drho = [drho; dphi];                 % Complete relative dynamics vector field
    end
end

%% Auxiliary functions 
% Geometrical regularized integration 
function [drho] = geometrical_model(mu, omega, s_t, s_r)
    % Constants of the system 
    mu_r(1) = 1-mu;               % Reduced gravitational parameter of the first primary 
    mu_r(2) = mu;                 % Reduced gravitational parameter of the second primary 

    % Synodic position of the primaries 
    R(:,1) = [-mu; 0; 0];         % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];        % Synodic position of the second primary
    
    % State variables 
    r_r = s_r(1:3);               % Synodic relative position 
    v_r = s_r(4:6);               % Synodic relative velocity
    r_t = s_t(1:3);               % Synodic position of the target
        
    % Encke acceleration method
    gamma = [2*v_r(2)/omega+r_r(1)/omega^2; -2*v_r(1)/omega+r_r(2)/omega^2; 0];
    for i = 1:length(mu_r)
        q = -dot(2*(r_t-R(:,i))+r_r,r_r)/norm(r_t+r_r-R(:,i))^2;
        f = q*(3*(1+q)+q^2)/(1+(1+q)^(3/2));
        gamma = gamma - (mu_r(i)/norm(r_t-R(:,i))^3)*(f*(r_t-R(:,i))+(1+f)*r_r)/omega^2;
    end
    
    % Equations of motion 
    drho = [v_r; 
            gamma];
end