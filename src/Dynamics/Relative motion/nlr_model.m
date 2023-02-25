%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 10/03/21
% File: fulrel_motion.m 
% Issue: 0 
% Validated: 

%% Full relative motion in the CR3BP %%
% This function contains the full nonlinear relative motion between two punctual particles in the 
% CR3BP. It accounts for the two masses moving in the normalized, non dimensional synodic frame 
% define by the two primaries, which are assumed to be in the same  plane and in circular orbits

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - scalar direction (in binary format, 1 or -1), indicating the
%           time integration direction: 1 for forward integration, -1 for
%           backward integration
%         - boolean flagVar, true for dyanmics and STM integration, 
%           false for only dynamical integration
%         - boolean relFlagVar, true for relative dyanmics and STM integration, 
%           false for only dynamical integration
%         - string method_ID, identifying which integration method to use
%         - scalar t, a reference epoch
%         - vector s, containing in an Nx1 array the phase space vector,
%           possibly augmented with the state transition matrix at time t
%         - cell array varargin, to include GNC requirements on the motion
%           of the spacecraft

% Outputs: - vector ds, the differential vector field, which will include
%            the phase space trajectory

% Methods: 

% New versions: 

function [ds] = nlr_model(mu, direction, flagVar, relFlagVar, method_ID, t, s, varargin)
    % Constants 
    m = 6;                                                        % Individual phase space dimension
    
    % State variables 
    if (flagVar)
        s_t = s(1:m+m^2);                                         % State of the target
        s_r = s(m+m^2+1:end);                                     % State of the relative particle
    else
        s_t = s(1:m);                                             % State of the target
        s_r = s(7:end);                                           % State of the relative particle  
    end

    % Equations of motion of the target
    if (~isempty(varargin))
        GNC = varargin{1};                                        % GNC handling structure
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
    s_t = s_t(1:m);                                               % Eliminate the variational state

    
    % Equations of motion of the relative state 
    switch (method_ID)
        % Deterministic models
        case 'Encke'
            drho = Encke_method(mu, s_t, s_r, varargin);          % Relative motion equations
        case 'Full nonlinear'
            drho = full_model(mu, s_t, s_r);                      % Relative motion equations
        case 'Linear' 
            drho = rlm_model(mu, s_t, s_r);                       % Relative motion equations
        case 'Second order' 
            drho = second_order_model(mu, s_t, s_r);              % Relative motion equations 
        case 'Third order'
            drho = third_order_model(mu, s_t, s_r);               % Relative motion equations

        % Uncertainty models 
        case 'Two body'
            drho = two_body(mu, s_t, s_r);
        otherwise
            error('No valid model was chosen');
    end

    % Augmented integral state dynamics
    if (mod(length(s_r), 9) == 0)
        m = 9;                                      % Augmented phase space dimension
        drho = [drho; s_r(7:9)];                    % Integrate the relative position for the PID controller
    end
    
    % GNC handler 
    if (~isempty(GNC))
        % GNC scheme
        [~, ~, u] = GNC_handler(GNC, s_t(1:6).', s_r.', t);         % Compute the control law
        drho(4:6) = drho(4:6)+u;                                         % Add the control vector       
    end
    
    % Relative variational equations
    if (relFlagVar)
        dphi = var_equations(mu, t, s);      % Linear variational equations
        drho = [drho; dphi];                 % Complete relative dynamics vector field
    end
    
    % Vector field 
    ds = [ds_t; drho];
end

%% Auxiliary functions 
%Full nonlinear relative motion equations
function [drho] = full_model(mu, s_t, s_r)
    % Constants of the system 
    mup(1) = 1-mu;              % Reduced gravitational parameter of the first primary 
    mup(2) = mu;                % Reduced gravitational parameter of the second primary 
    
    % State variables 
    r_t = s_t(1:3);             % Synodic position of the target
    r_r = s_r(1:3);             % Synodic relative position 
    v_r = s_r(4:6);             % Synodic relative velocity
    
    x = r_r(1);                 % Synodic relative x coordinate
    y = r_r(2);                 % Synodic relative y coordinate
    
    % Synodic position of the primaries 
    R(:,1) = [-mu; 0; 0];       % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];      % Synodic position of the second primary
    Rr(:,1) = r_t-R(:,1);       % Synodic relative position of the target to the first primary
    Rr(:,2) = r_t-R(:,2);       % Synodic relative position of the target to the second primary
    
    % Relative acceleration 
    gamma = [2*v_r(2)+x; -2*v_r(1)+y; 0];                                            % Synodic acceleration
    F(:,1) = mup(1)*(Rr(:,1)/norm(Rr(:,1))^3-(r_r+Rr(:,1))/norm(r_r+Rr(:,1))^3);     % Gravitational force of the first primary
    F(:,2) = mup(2)*(Rr(:,2)/norm(Rr(:,2))^3-(r_r+Rr(:,2))/norm(r_r+Rr(:,2))^3);     % Gravitational force of the second primary
    gamma = gamma + sum(F,2);                                                        % Total synodic acceleration
    
    %Equations of motion 
    drho = [v_r; 
            gamma];
end

% Relative motion equations linearized with respect to the target
function [drho] = rlm_model(mu, s_t, s_r)
    % Constants of the system 
    mup(1) = 1-mu;                          % Reduced gravitational parameter of the first primary 
    mup(2) = mu;                            % Reduced gravitational parameter of the second primary 
    
    % State variables 
    r_t = s_t(1:3);                         % Synodic position of the target
    rho = s_r(1:3);                         % Relative position vecto

    % Synodic position of the primaries 
    R(:,1) = [-mu; 0; 0];                   % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];                  % Synodic position of the second primary
    
    % Relative position between the primaries and the target 
    Ur(:,1) = R(:,1)-r_t;                   % Position of the target with respect to the first primary
    ur(:,1) = Ur(:,1)/norm(Ur(:,1));        % Unit vector of the relative position of the target with respect to the first primary
    Ur(:,2) = R(:,2)-r_t;                   % Position of the target with respect to the first primary
    ur(:,2) = Ur(:,2)/norm(Ur(:,2));        % Unit vector of the relative position of the target with respect to the second primary
    
    % Relative acceleration vector field
    O = zeros(3);                           % 3 by 3 null matrix
    I = eye(3);                             % 3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];         % Hat map dyadic of the angular velocity for the synodice reference frame
    
    % Relative acceleration (linear order term)
    Sigma = zeros(3);
    for i = 1:length(mup)
        % Relative Legendre coefficients          
        c2 = mup(i) / norm(Ur(:,i))^3;      % Third order relative Legendre coefficient
        
        % Compute the acceleration
        Sigma = Sigma + c2 * ( -I + 3 * ur(:,i)*ur(:,i).' );
    end
    
    % State matrix 
    A = [O I; Sigma 2*Omega];
    
    % Equations of motion 
    drho = A*s_r;
end

% Second order relative motion equations 
function [drho] = second_order_model(mu, s_t, s_r)  
    % Constants of the system 
    mup(1) = 1-mu;              % Reduced gravitational parameter of the first primary 
    mup(2) = mu;                % Reduced gravitational parameter of the second primary

    % State variables 
    r_t = s_t(1:3);             % Position vector of the target
    rho = s_r(1:3);             % Relative position vector

    % Synodic position of the primaries 
    R(:,1) = [-mu; 0; 0];       % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];      % Synodic position of the second primary
    Rr(:,1) = R(:,1)-r_t;       % Synodic relative position of the target to the first primary
    Rr(:,2) = R(:,2)-r_t;       % Synodic relative position of the target to the second primary
            
    % Relative acceleration (linear term)
    drho = rlm_model(mu, s_t, s_r);
    
    % Relative acceleration (third order term)
    Sigma = dot(rho,rho)*eye(3)-rho*rho.';
    gamma = zeros(6,1);
    for i = 1:length(mup)
        % Relative Legendre coefficients          
        c3 = mup(i) / norm(Rr(:,i))^4;      % Third order relative Legendre coefficient
        
        % Compute the acceleration
        cos_theta = dot(rho, Rr(:,i))/(norm(rho)*norm(Rr(:,i)));
        gamma(4:6) = c3/2 * (3 * norm(rho) * (-3*cos_theta) * rho + (15*cos_theta^2-3) * Sigma.' * Rr(:,i)/norm(R(:,i)) );
        drho = drho + gamma;
    end
end

% Third order relative motion equations 
function [drho] = third_order_model(mu, s_t, s_r)  
    % Constants of the system 
    mup(1) = 1-mu;              % Reduced gravitational parameter of the first primary 
    mup(2) = mu;                % Reduced gravitational parameter of the second primary

    % State variables 
    r_t = s_t(1:3);             % Position vector of the target
    rho = s_r(1:3);             % Relative position vector

    % Synodic position of the primaries 
    R(:,1) = [-mu; 0; 0];       % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];      % Synodic position of the second primary
    Rr(:,1) = R(:,1)-r_t;       % Synodic relative position of the target to the first primary
    Rr(:,2) = R(:,2)-r_t;       % Synodic relative position of the target to the second primary
            
    % Relative acceleration (linear and second order term)
    drho = second_order_model(mu, s_t, s_r);
    
    % Relative acceleration (third order term)
    Sigma = dot(rho,rho)*eye(3)-rho*rho.';
    gamma = zeros(6,1);
    for i = 1:length(mup)
        % Relative Legendre coefficients          
        c4 = mup(i) / norm(Rr(:,i))^5;      % Third order relative Legendre coefficient

        % Compute the acceleration
        cos_theta = dot(rho, Rr(:,i))/(norm(rho)*norm(Rr(:,i)));
        gamma(4:6) = c4/8 * (4 * norm(rho)^2 * (-30*cos_theta^2+3) * rho + norm(rho) * (140*cos_theta^3-60*cos_theta^2) * Sigma.' * Rr(:,i)/norm(R(:,i)) );
        drho = drho + gamma;
    end
end

% Full nonlinear relative motion equations via Encke's method
function [drho] = Encke_method(mu, s_t, s_r, varargin)
    % Constants of the system 
    mu_r(1) = 1-mu;               % Reduced gravitational parameter of the first primary 
    mu_r(2) = mu;                 % Reduced gravitational parameter of the second primary 
    
    % State variables 
    r_t = s_t(1:3);               % Synodic position of the target
    r_r = s_r(1:3);               % Synodic relative position 
    v_r = s_r(4:6);               % Synodic relative velocity 
    
    % Synodic position of the primaries 
    R(:,1) = [-mu; 0; 0];         % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];        % Synodic position of the second primary
    
    % Encke acceleration method
    gamma = [2*v_r(2)+r_r(1); -2*v_r(1)+r_r(2); 0];
    for i = 1:length(mu_r)
        q = -dot(2*(r_t-R(:,i))+r_r,r_r)/norm(r_t+r_r-R(:,i))^2;
        f = q*(3*(1+q)+q^2)/(1+(1+q)^(3/2));
        gamma = gamma - (mu_r(i)/norm(r_t-R(:,i))^3)*(f*(r_t-R(:,i))+(1+f)*r_r);
    end
    
    % Equations of motion 
    drho = [v_r; 
            gamma];
end

% Two body solution around the second primary as a first order approximation of the problem, with a numerical Encke's method
function [drho] = two_body(mu, s_r, s_t)
    % State variables 
    r_t = s_t(1:3);               % Synodic position of the target
    r_r = s_r(1:3);               % Synodic relative position 
    v_r = s_r(4:6);               % Synodic relative velocity 
        
    % Encke acceleration method
    gamma = [2*v_r(2)+r_r(1); -2*v_r(1)+r_r(2); 0];
    q = -dot(2*r_t+r_r,r_r)/norm(r_t+r_r)^2;
    f = q*(3*(1+q)+q^2)/(1+(1+q)^(3/2));
    gamma = gamma - (mu/norm(r_t)^3)*(f*(r_t)+(1+f)*r_r);
    
    %Equations of motion 
    drho = [v_r; 
            gamma];
end