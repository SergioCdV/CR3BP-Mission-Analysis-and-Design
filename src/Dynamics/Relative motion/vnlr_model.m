%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 03/10/22
% File: vnlr_model.m 
% Issue: 0 
% Validated: 

%% Full vectorized relative motion in the CR3BP %%
% This function contains the full nonlinear vectorized relative motion between two punctual particles in the 
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

% New versions: vectorize all models and options

% Full nonlinear relative motion equations via Encke's method
function [drho] = vnlr_model(mu, u, t, s)
    % Constants of the system 
    n = 6;                        % Phase space dimension
    mu_r(1) = 1-mu;               % Reduced gravitational parameter of the first primary 
    mu_r(2) = mu;                 % Reduced gravitational parameter of the second primary 

    s_t = s(1:n,:);               % Target state
    s_r = s(n+1:end,:);           % Chaser relative state
        
    % Synodic position of the primaries 
    R(:,1) = [-mu; 0; 0];                                   % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];                                  % Synodic position of the second primary

    % Target differential vectorfield
    Rr(1:3,:) = s_t(1:3,:)-R(:,1);                          % Relative position vector to the first primary
    Rr(4:6,:) = s_t(1:3,:)-R(:,2);                          % Relative position vector to the secondary primary
    
    % Compute the time flow of the system (depending on the time direction)
    omega = [2*s_t(5,:)+s_t(1,:); -2*s_t(4,:)+s_t(2,:); zeros(1,size(s_t,2))];
    dst = [s_t(4:6,:); omega-mu_r(1)./sqrt(dot(Rr(1:3,:),Rr(1:3,:),1)).^3.*Rr(1:3,:)-mu_r(2)./sqrt(dot(Rr(4:6,:),Rr(4:6,:),1)).^3.*Rr(4:6,:)];       % Time flow of the system

    % Relative dynamics constants
    R = -Rr;
    Rr(1:3,:) = s_r(1:3,:)-R(1:3,:);             % Relative position to the first primary
    Rr(4:6,:) = s_r(1:3,:)-R(4:6,:);             % Relative position to the first primary

    % First primary gravitational acceleration
    q = -dot(s_r(1:3,:)-2*R(1:3,:),s_r(1:3,:),1)./sqrt(dot(Rr(1:3,:),Rr(1:3,:),1)).^2;
    f = q.*(3*(1+q)+q.^2)./(1+(1+q).^(3/2));
    gamma(1:3,:) = -mu_r(1)./sqrt(dot(R(1:3,:),R(1:3,:),1)).^3.*((1+f).*s_r(1:3,:)-f.*R(1:3,:));

    % Second primary gravitational acceleration
    q = -dot(s_r(1:3,:)-2*R(4:6,:),s_r(1:3,:),1)./sqrt(dot(Rr(4:6,:),Rr(4:6,:),1)).^2;
    f = q.*(3*(1+q)+q.^2)./(1+(1+q).^(3/2));
    gamma(4:6,:) = -mu_r(2)./sqrt(dot(R(4:6,:),R(4:6,:),1)).^3.*((1+f).*s_r(1:3,:)-f.*R(4:6,:));

    % Gyroscopic acceleration
    omega = [2*s_r(5,:)+s_r(1,:); -2*s_r(4,:)+s_r(2,:); zeros(1,size(s_r,2))];
    
    % Final vectorfield
    dv = s_r(4:6,:);                              % Inertial velocity field
    f = gamma(1:3,:)+gamma(4:6,:)+omega;          % Forces field
    drho = [dst;
            dv; 
            f+u];                                 % Inertial acceleration field 

    if (mod(size(s,1)-n,n+n/2) == 0)
        n = n + n/2;                              % Augmented phase space dimension
        drho = [drho; s_r(1:3,:)];                % Augmented dynamics
    end

    if (size(s_r,1) == n+n^2)
        dphi = var_equations(mu, t, s);
        drho = [drho; dphi];
    end
end

function [ds] = to_be_model(mu, direction, flagVar, relFlagVar, method_ID, t, s, varargin)
    % Constants 
    m = 6;                                                        % Individual phase space dimension
    
    % State variables 
    if (flagVar)
        s_t = s(1:m+m^2);                                         % State of the target
        s_r = s(m+m^2+1:end);                                     % State of the relative particle
    else
        s_t = s(1:6);                                             % State of the target
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
    s_t = s_t(1:6);                                               % Eliminate the variational state
    
    % Equations of motion of the relative state 
    switch (method_ID)
        % Deterministic models
        case 'Encke'
            drho = Encke_method(mu, s_t, s_r, varargin);          % Relative motion equations
        case 'Full nonlinear'
            drho = full_model(mu, s_t, s_r);                      % Relative motion equations
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
            case 'HDRE'
                drho = [drho; s_r(7:9)];
            case 'SMC'
            case 'RFSK'
            case 'MFSK'
            case 'PFSK'
            case 'TAHO'
            case 'iLQR'
            otherwise
                error('No valid control algorithm was selected');
        end
        
        % GNC scheme
        [~, ~, u] = GNC_handler(GNC, s_t(1:6).', s_r.', t);         % Compute the control law
        drho(4:6) = drho(4:6)+u;                                    % Add the control vector       
    end
    
    % Relative variational equations
    if (relFlagVar)
        Phi = reshape(s_r(7:end), [m m]);               % State Transition Matrix
        J = rel_jacobian(mu, [s_t; s_r(1:3)]);          % Relative Jacobian matrix
        dphi = J*Phi;                                   % New State Transition Matrix
        dphi = reshape(dphi, [m^2 1]);                  % First linear variational vector field
        drho = [drho; dphi];                            % Complete relative dynamics vector field
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
    Rr(:,2) = r_t-R(:,1);       % Synodic relative position of the target to the second primary
    
    % Relative acceleration 
    gamma = [2*v_r(2)+y; -2*v_r(1)+x; 0];                                            % Synodic acceleration
    F(:,1) = mup(1)*(Rr(:,1)/norm(Rr(:,1))^3-(r_r+Rr(:,1))/norm(r_r+Rr(:,1))^3);     % Gravitational force of the first primary
    F(:,2) = mup(2)*(Rr(:,2)/norm(Rr(:,2))^3-(r_r+Rr(:,2))/norm(r_r+Rr(:,2))^3);     % Gravitational force of the second primary
    gamma = gamma + sum(F,2);                                                        % Total synodic acceleration
    
    %Equations of motion 
    drho = [v_r; 
            gamma];
end

% Second order relative motion equations 
function [drho] = second_order_model(mu, s_t, s_r)  
    % State variables 
    r_t = s_t(1:3);                             % Position vector of the target
    x = s_r(1);                                 % Synodic x coordinate of the relative position
    y = s_r(2);                                 % Synodic y coordinate of the relative position
    z = s_r(3);                                 % Synodic z coordinate of the relative position
        
    % Relative Legendre coefficients          
    cn = relegendre_coefficients(mu, r_t, 3);   % Relative Legendre coefficients 
    c2 = cn(2);                                 % First order relative Legendre coefficient
    c3 = cn(3);                                 % Second order relative Legendre coefficient
    
    % Relative acceleration
    O = zeros(3,3);                             % 3 by 3 null matrix
    I = eye(3);                                 % 3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];             % Hat map dyadic of the angular velocity for the synodice reference frame
    
    % Gravity acceleration
    Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];                           % Linear term
    Sigma3 = c3*[0; 0; 0; (3/2)*(2*x^2-y^2-z^2); -3*x*y; -3*x*z];      % Second order term
    
    % Equations of motion 
    drho = [O I; Sigma 2*Omega]*s_r + Sigma3;
end

% Third order relative motion equations 
function [drho] = third_order_model(mu, s_t, s_r)  
    % State variables 
    r_t = s_t(1:3);                             % Position vector of the target
    x = s_r(1);                                 % Synodic x coordinate of the relative position
    y = s_r(2);                                 % Synodic y coordinate of the relative position
    z = s_r(3);                                 % Synodic z coordinate of the relative position
        
    % Relative Legendre coefficients          
    cn = relegendre_coefficients(mu, r_t, 4);   % Relative Legendre coefficients 
    c2 = cn(2);                                 % First order relative Legendre coefficient
    c3 = cn(3);                                 % Second order relative Legendre coefficient
    c4 = cn(4);                                 % Third order relative Legendre coefficient
    
    % Relative acceleration (non inertial)
    O = zeros(3,3);                             % 3 by 3 null matrix
    I = eye(3);                                 % 3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];             % Hat map dyadic of the angular velocity for the synodice reference frame
    
    % Gravity acceleration
    Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];                                    % Linear term
    Sigma3 = c3*[0; 0; 0; (3/2)*(2*x^2-y^2-z^2); -3*x*y; -3*x*z];               % Second order term
    Sigma4 = c4*[0; 0; 0; 2*x*(2*x^2-3*y^2-3*z^2); ...
                 -(3/2)*y*(4*x^2-y^2-z^2); -(3/2)*z*(4*x^2-y^2-z^2)];           % Third order term
    
    %Equations of motion 
    drho = [O I; Sigma 2*Omega]*s_r + Sigma3 + Sigma4;
end