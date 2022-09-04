%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/20
% File: rel_motion.m 
% Issue: 0 
% Validated: 

%% Linear relative motion model in the CR3BP %%
% This function contains several linear relative motion models between two punctual particles in the 
% CR3BP. It accounts for the two masses moving in the normalized, non dimensional synodic frame define 
% by the two primaries, which are assumed to be in the same plane and in circular orbits

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - vector cn, with the Legendre coefficients of the relative
%           motion
%         - scalar direction (in binary format, 1 or -1), indicating the
%           time integration direction: 1 for forward integration, -1 for
%           backward integration
%         - boolean flagVar, true for dyanmics and STM integration, 
%           false for only dynamical integration
%         - string model, indicating the linear model to be used
%         - scalar t, a reference epoch
%         - vector s, containing in an Nx1 array the phase space vector,
%           possibly augmented with the state transition matrix at time t
%         - cell array varargin, to include GNC requirements on the motion
%           of the spacecraft

% Outputs: - vector dr, the differential vector field, which will include
%            the phase space trajectory

% Methods: different linear models of relative motion

% New versions: include the first variations of the vector field

function [ds] = lr_model(mu, cn, direction, flagVar, model, t, s, varargin)
    % State variables 
    s_t = s(1:6);                                % State of the target
    s_r = s(7:12);                               % State of the chaser
    
    % Equations of motion of the target
    if (~isempty(varargin))
        GNC = varargin{1};                       % GNC handling structure
        if (iscell(GNC))
            GNC = GNC{1};
        end
        % Target equations of motion
        if (isfield(GNC, 'Target'))
            ds_t = cr3bp_equations(mu, direction, flagVar, t, s_t, GNC.Target); 
        else
            ds_t = cr3bp_equations(mu, direction, flagVar, t, s_t); 
        end
    else
        GNC = [];                                                 % Empty GNC structure
        ds_t = cr3bp_equations(mu, direction, flagVar, t, s_t);   % Target equations of motion
    end
    s_t = s_t(1:6);                                               % Eliminate the variational state                                               

    % Equations of motion of the relative state 
    switch (model)
        case 'RLM'
            drho = rlm_model(mu, s_t, s_r);                       % Relative motion equations
        case 'SLLM'
            drho = sllm_model(cn, s_r);                           % Relative motion equations
        case 'ULLM'
            drho = ullm_model(mu, s_t, s_r);                      % Relative motion equations
        otherwise
            error('No valid linear model was chosen');
    end
    
    % GNC handler 
    if (~isempty(GNC))        
        % Integrate the relative position for the PID controller
        switch (GNC.Algorithms.Control)
            case 'TISS'
                error('No valid control algorithm was selected for integration purposes')
            case 'MISS'
                error('No valid control algorithm was selected for integration purposes')
            case 'TITA'
                error('No valid control algorithm was selected for integration purposes')
            case 'MPC'
                error('No valid control algorithm was selected for integration purposes')
            case 'APF'
                error('No valid control algorithm was selected for integration purposes')
            case 'SDRE'
                drho = [drho; s_r(7:9)];
            case 'LQR'
                drho = [drho; s_r(7:9)];
            case 'SMC'
                
            otherwise
                error('No valid control algorithm was selected');
        end
        
        % GNC scheme
        [~, ~, u] = GNC_handler(GNC, s_t(1:6).', s_r.', t);         % Compute the control law
        drho(4:6) = drho(4:6)+u;                                    % Add the control vector       
    end
        
    % Vector field 
    ds = [ds_t; drho];
end

%% Auxiliary functions 
% Relative motion equations linearized with respect to the target
function [drho] = rlm_model(mu, s_t, s_r)
    % Constants of the system 
    mup(1) = 1-mu;                          % Reduced gravitational parameter of the first primary 
    mup(2) = mu;                            % Reduced gravitational parameter of the second primary 
    
    % State variables 
    r_t = s_t(1:3);                         % Synodic position of the target
    
    % Synodic position of the primaries 
    R(:,1) = [-mu; 0; 0];                   % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];                  % Synodic position of the second primary
    
    % Relative position between the primaries and the target 
    Ur(:,1) = r_t-R(:,1);                   % Position of the target with respect to the first primary
    ur(:,1) = Ur(:,1)/norm(Ur(:,1));        % Unit vector of the relative position of the target with respect to the first primary
    Ur(:,2) = r_t-R(:,2);                   % Position of the target with respect to the first primary
    ur(:,2) = Ur(:,2)/norm(Ur(:,2));        % Unit vector of the relative position of the target with respect to the second primary
    
    % Relative acceleration vector field
    O = zeros(3);                           % 3 by 3 null matrix
    I = eye(3);                             % 3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];         % Hat map dyadic of the angular velocity for the synodice reference frame
    
    % Gravity acceleration
    Sigma = -((mup(1)/norm(Ur(:,1))^3)+(mup(2)/norm(Ur(:,2)))^3)*eye(3) ...
            +3*((mup(1)/norm(Ur(:,1))^3)*(ur(:,1)*ur(:,1).')+(mup(2)/norm(Ur(:,2))^3)*(ur(:,2)*ur(:,2).'));
    
    % State matrix 
    A = [O I; Sigma-Omega*Omega 2*Omega];
    
    % Equations of motion 
    drho = A*s_r;
end

% Relative motion equations linearized with respect to the absolute libration point
function [drho] = sllm_model(cn, s_r)    
    % Legendre coefficient c2           
    c2 = cn(2);                                 % Second Legendre coefficient
    
    % Relative acceleration (non inertial)
    O = zeros(3);                               % 3 by 3 null matrix
    I = eye(3);                                 % 3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];             % Hat map dyadic of the angular velocity for the synodice reference frame
    Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];    % Gravity acceleration
    
    %Equations of motion 
    drho = [O I; Sigma 2*Omega]*s_r;
end

% Relative motion equations linearized with respect to the relative libration point
function [drho] = ullm_model(mu, s_t, s_r)  
    % State variables 
    r_t = s_t(1:3);                             % Position vector of the target
        
    % Relative Legendre coefficient c2           
    cn = relegendre_coefficients(mu, r_t, 2);   % Relative Legendre coefficients 
    c2 = cn(2);                                 % First order relative Legendre coefficient
    
    % Relative acceleration (non inertial)
    O = zeros(3);                               % 3 by 3 null matrix
    I = eye(3);                                 % 3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];             % Hat map dyadic of the angular velocity for the synodice reference frame
    Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];    % Gravity acceleration
    
    % Equations of motion 
    drho = [O I; Sigma 2*Omega]*s_r;
end
