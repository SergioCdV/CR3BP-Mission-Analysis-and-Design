%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 10/03/21
% File: fulrel_motion.m 
% Issue: 0 
% Validated: 

%% Full relative motion in the CR3BP %%
% This function contains the full nonlinear relative motion between two punctual particles in the 
% circular restricted three body problem. It accounts for the two masses moving in the normalized, 
% non dimensional synodic frame define by the two primaries, which are assumed to be in the same 
% plane and in circular orbits. 

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - scalar direction (in binary format, 1 or -1), indicating the
%           time integration direction: 1 for forward integration, -1 for
%           backward integration.
%         - boolean flagVar, true for dyanmics and STM integration, 
%           false for only dynamical integration.
%         - string method_ID, identifying which integration method to use.
%         - scalar t, a reference epoch. 
%         - vector s, containing in an Nx1 array the phase space vector,
%           possibly augmented with the state transition matrix at time t. 

% Outputs: - vector ds, the differential vector field, which will include
%            the phase space trajectory.

% Methods: . 

% New versions: 

function [ds] = nlr_model(mu, direction, flagVar, relFlagVar, method_ID, t, s, varargin)
    %Constants 
    m = 6;                                                        %Individual phase space dimension
    
    %State variables 
    if (flagVar)
        s_t = s(1:m+m^2);                                         %State of the target
        s_r = s(m+m^2+1:end);                                     %State of the relative particle
    else
        s_t = s(1:6);                                             %State of the target
        s_r = s(7:end);                                           %State of the relative particle  
    end

    %Equations of motion of the target
    if (~isempty(varargin))
        GNC = varargin{1};                                        %GNC handling structure
        if (iscell(GNC))
            GNC = GNC{1};
        end
        %Target equations of motion
        if (isfield(GNC, 'Target'))
            ds_t = cr3bp_equations(mu, direction, flagVar, t, s_t, GNC.Target); 
        else
            ds_t = cr3bp_equations(mu, direction, flagVar, t, s_t); 
        end
    else
        GNC = [];                                                 %Empty GNC structure
        ds_t = cr3bp_equations(mu, direction, flagVar, t, s_t);   %Target equations of motion
    end
    s_t = s_t(1:6);                                               %Eliminate the variational state
    
    %Equations of motion of the relative state 
    switch (method_ID)
        case 'Encke'
            drho = Encke_method(mu, s_t, s_r, varargin);          %Relative motion equations
        case 'Full nonlinear'
            drho = full_model(mu, s_t, s_r);                      %Relative motion equations
        case 'Second order' 
            drho = second_order_model(mu, s_t, s_r);              %Relative motion equations 
        case 'Third order'
            drho = third_order_model(mu, s_t, s_r);               %Relative motion equations
        otherwise
            error('No valid model was chosen');
    end
    
    %GNC handler 
    if (~isempty(GNC))        
        %Integrate the relative position for the PID controller
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
        
        %GNC scheme
        [~, ~, u] = GNC_handler(GNC, s_t(1:6).', s_r.', t);         %Compute the control law
        drho(4:6) = drho(4:6)+u;                                    %Add the control vector       
    end
    
    %Relative variational equations
    if (relFlagVar)
        Phi = reshape(s_r(7:end), [m m]);               %State transition matrix
        J = rel_jacobian(mu, [s_t; s_r(1:3)]);          %Relative Jacobian matrix
        dphi = J*Phi;                                   %New state transition matrix
        dphi = reshape(dphi, [m^2 1]);                  %Differential vector on the state transition matrix
        drho = [drho; dphi];                            %Complete relative dynamics vector field
    end
    
    %Vector field 
    ds = [ds_t; drho];
end

%% Auxiliary functions 
%Full nonlinear relative motion equations
function [drho] = full_model(mu, s_t, s_r)
    %Constants of the system 
    mu1 = 1-mu;             %Reduced gravitational parameter of the first primary 
    mu2 = mu;               %Reduced gravitational parameter of the second primary 
    
    %State variables 
    r_t = s_t(1:3);         %Synodic position of the target
    r_r = s_r(1:3);         %Synodic relative position 
    v_r = s_r(4:6);         %Synodic relative velocity
    
    x = r_r(1);             %Synodic relative x coordinate
    y = r_r(2);             %Synodic relative y coordinate
    
    %Synodic position of the primaries 
    R1 = [-mu; 0; 0];       %Synodic position of the first primary
    R2 = [1-mu; 0; 0];      %Synodic position of the second primary
    
    %Relative acceleration 
    gamma = [2*v_r(2)+y; -2*v_r(1)+x; 0];                                   %Synodic acceleration
    F1 = mu1*((r_t-R1)/norm(r_t-R1)^3-(r_t-R1+r_r)/norm(r_t-R1+r_r)^3);     %Gravitational force of the first primary
    F2 = mu2*((r_t-R2)/norm(r_t-R2)^3-(r_t-R2+r_r)/norm(r_t-R2+r_r)^3);     %Gravitational force of the second primary
    gamma = gamma + F1 + F2;                                                %Total synodic acceleration
    
    %Equations of motion 
    drho = [v_r; 
            gamma];
end

%Second order relative motion equations 
function [drho] = second_order_model(mu, s_t, s_r)  
    %State variables 
    r_t = s_t(1:3);                             %Position vector of the target
    x = s_r(1);                                 %Synodic x coordinate of the relative position
    y = s_r(2);                                 %Synodic y coordinate of the relative position
    z = s_r(3);                                 %Synodic z coordinate of the relative position
        
    %Relative Legendre coefficients          
    cn = relegendre_coefficients(mu, r_t, 3);   %Relative Legendre coefficients 
    c2 = cn(2);                                 %First order relative Legendre coefficient
    c3 = cn(3);                                 %Second order relative Legendre coefficient
    
    %Relative acceleration (non inertial)
    O = zeros(3,3);                             %3 by 3 null matrix
    I = eye(3);                                 %3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];             %Hat map dyadic of the angular velocity for the synodice reference frame
    
    %Gravity acceleration
    Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];                           %Linear term
    Sigma3 = c3*[0; 0; 0; (3/2)*(2*x^2-y^2-z^2); -3*x*y; -3*x*z];      %Second order term
    
    %Equations of motion 
    drho = [O I; Sigma -2*Omega]*s_r + Sigma3;
end

%Third order relative motion equations 
function [drho] = third_order_model(mu, s_t, s_r)  
    %State variables 
    r_t = s_t(1:3);                             %Position vector of the target
    x = s_r(1);                                 %Synodic x coordinate of the relative position
    y = s_r(2);                                 %Synodic y coordinate of the relative position
    z = s_r(3);                                 %Synodic z coordinate of the relative position
        
    %Relative Legendre coefficients          
    cn = relegendre_coefficients(mu, r_t, 4);   %Relative Legendre coefficients 
    c2 = cn(2);                                 %First order relative Legendre coefficient
    c3 = cn(3);                                 %Second order relative Legendre coefficient
    c4 = cn(4);                                 %Third order relative Legendre coefficient
    
    %Relative acceleration (non inertial)
    O = zeros(3,3);                             %3 by 3 null matrix
    I = eye(3);                                 %3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];             %Hat map dyadic of the angular velocity for the synodice reference frame
    
    %Gravity acceleration
    Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];                                    %Linear term
    Sigma3 = c3*[0; 0; 0; (3/2)*(2*x^2-y^2-z^2); -3*x*y; -3*x*z];               %Second order term
    Sigma4 = c4*[0; 0; 0; 2*x*(2*x^2-3*y^2-3*z^2); ...
                 -(3/2)*y*(4*x^2-y^2-z^2); -(3/2)*z*(4*x^2-y^2-z^2)];           %Third order term
    
    %Equations of motion 
    drho = [O I; Sigma -2*Omega]*s_r + Sigma3 + Sigma4;
end

%Full nonlinear relative motion equations via Encke's method
function [drho] = Encke_method(mu, s_t, s_r, varargin)
    %Constants of the system 
    mu_r(1) = 1-mu;               %Reduced gravitational parameter of the first primary 
    mu_r(2) = mu;                 %Reduced gravitational parameter of the second primary 
    
    %State variables 
    r_t = s_t(1:3);               %Synodic position of the target
    r_r = s_r(1:3);               %Synodic relative position 
    v_r = s_r(4:6);               %Synodic relative velocity 
    
    %Synodic position of the primaries 
    R(1:3,1) = [-mu; 0; 0];       %Synodic position of the first primary
    R(1:3,2) = [1-mu; 0; 0];      %Synodic position of the second primary
    
    %Encke acceleration method
    gamma = [2*v_r(2)+r_r(1); -2*v_r(1)+r_r(2); 0];
    for i = 1:length(mu_r)
        q = -dot(2*(r_t-R(:,i))+r_r,r_r)/norm(r_t+r_r-R(:,i))^2;
        f = q*(3*(1+q)+q^2)/(1+(1+q)^(3/2));
        gamma = gamma - (mu_r(i)/norm(r_t-R(:,i))^3)*(f*(r_t-R(:,i))+(1+f)*r_r);
    end
    
    %Equations of motion 
    drho = [v_r; 
            gamma];
end