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

function [ds] = nlr_model(mu, direction, flagVar, method_ID, t, s, varargin)
    %Constants 
    m = 6;                                                        %Individual phase space dimension
    
    %State variables 
    if (flagVar)
        method_ID = 'Encke V';
        s_t = s(1:m+m^2);                                         %State of the target
    else
        s_t = s(1:6);                                             %State of the target
        s_r = s(7:12);                                            %State of the relative particle
    end
    
    %Equations of motion of the target
    ds_t = cr3bp_equations(mu, direction, flagVar, t, s_t);       %Target equations of motion
    
    %Equations of motion of the relative state 
    switch (method_ID)
        case 'Encke'
            drho = Encke_method(mu, s_t, s_r);                    %Relative motion equations
        case 'Encke V'
            drho = EnckeV_method(mu, flagVar, s);                 %Relative motion equations
        case 'Encke C'
            drho = EnckeC_method(mu, s, varargin);                %Relative motion equations
        case 'Encke LQR'
            drho = EnckeLQR_method(mu, s, varargin);              %Relative motion equations
        case 'Encke SDRE'    
            drho = EnckeSDRE_method(mu, s, varargin);             %Relative motion equations
        case 'Full nonlinear'
            drho = full_model(mu, s_t, s_r);                      %Relative motion equations
        case 'Second order'
            drho = so_model(mu, s_t, s_r);                        %Relative motion equations 
        case 'Third order'
            drho = th_model(mu, s_t, s_r);                        %Relative motion equations
        otherwise
            drho = [];
            disp('No valid model was chosen');
    end
    
    %Vector field 
    ds = [ds_t; drho];
end

%% Auxiliary functions 
%Full nonlinear relative motion equations via Encke's method
function [drho] = Encke_method(mu, s_t, s_r)
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

%Full nonlinear relative motion equations and first variational system via Encke's method
function [drho] = EnckeV_method(mu, flagVar, s)
    %System parameters 
    m = 6;                                       %Phase space dimension
    
    %State variables 
    if (flagVar)
        s_t = s(1:m);                            %Target state 
        s_r = s(m+m^2+1:2*m+m^2);                %Relative state
        ds = reshape(s(2*m+m^2+1:end), [m m]);   %State transition matrix
    else
        s_t = s(1:6);                            %Target state 
        s_r = s(7:12);                           %Relative state
        ds = reshape(s(13:end), [m m]);          %State transition matrix
    end

    %Compute the integration of the relative motion equations 
    drho = Encke_method(mu, s_t, s_r);
    
    %Variational equations
    J = rel_jacobian(mu, s);            %Jacobian of the system
    ds = J*ds;                          %Variational equations 
    ds = reshape(ds, [m^2 1]);          %Variational vector field 
    
    %Final vector field 
    drho = [drho; ds]; 
end

%Full nonlinear relative motion equations with control vector
function [drho] = EnckeC_method(mu, s, varargin)
    %System parameters 
    m = 6;                                   %Phase space dimension
    
    %Control vector 
    aux = varargin{1};
    u = aux{1};
    
    %State variables
    s_t = s(1:m);                            %Target state 
    s_r = s(m+1:2*m);                        %Relative state

    %Compute the integration of the relative motion equations 
    drho = Encke_method(mu, s_t, s_r);                          %Natural vector field flow
    drho = drho + [0; 0; 0; u];                                 %Add control vector
end

%Full nonlinear relative motion equations with a continuous LQR PID controller
function [drho] = EnckeLQR_method(mu, s, varargin)
    %System parameters 
    m = 6;                                   %Phase space dimension
    
    %Control vector 
    aux = varargin{1};
    K = aux{1};
    
    %State variables
    s_t = s(1:m);                            %Target state 
    s_r = s(m+1:2*m);                        %Relative state
    int = s(end-2:end);                      %Integrator control vector
    
    %Compute the control law 
    u = -K*[int; s_r];

    %Compute the integration of the relative motion equations 
    drho = Encke_method(mu, s_t, s_r);                          %Natural vector field flow
    dint = s_r(1:3);                                            %Integrator field flow
    drho = drho + [0; 0; 0; u];                                 %Add control vector
    drho = [drho; dint];
end

%Full nonlinear relative motion equations with a continuous LQR PID controller
function [drho] = EnckeSDRE_method(mu, s, varargin)
    %System parameters 
    m = 6;                                  %Phase space dimension
    
    %Control vector 
    aux = varargin{1};
    model = aux{1};                         %Model to be used
    Ln = aux{2};                            %Libration point
    gamma = aux{3};                         %Relative distance of the libration point to the primary
    
    %State variables
    s_t = s(1:m);                           %Target state 
    r_t = s_t(1:3);                         %Target position
    s_r = s(m+1:2*m);                       %Relative state
    int = s(end-2:end);                     %Integrator control vector
    
    %Approximation 
    n = 6;                                  %Dimension of the state vector
    order = 2;                              %Order of the approximation 

    %Model coefficients 
    mup(1) = 1-mu;                          %Reduced gravitational parameter of the first primary 
    mup(2) = mu;                            %Reduced gravitational parameter of the second primary 
    R(:,1) = [-mu; 0; 0];                   %Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];                  %Synodic position of the second primary

    %Linear model matrices
    B = [zeros(n/2); zeros(n/2); eye(n/2)];         %Linear model input matrix 
    Omega = [0 2 0; -2 0 0; 0 0 0];                 %Coriolis dyadic
    
    %Cost function matrices 
    Q = diag(ones(1,n+3));                          %Cost weight on the state error
    M = eye(n/2);                                   %Cost weight on the spent energy

    %Select linear model 
    switch (model)
        case 'Fixed point'
            cn = legendre_coefficients(mu, Ln, gamma, order);     %Compute the relative Legendre coefficient c2 
            c2 = cn(2); 
            Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              %Linear model state matrix
        case 'Moving point' 
            rc = relegendre_coefficients(mu, r_t.', order);       %Compute the relative Legendre coefficient c2 
            c2 = rc(2); 
            Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              %Linear model state matrix
        case 'Target' 
            %Relative position between the primaries and the target 
            Ur1 = r_t-R(:,1);               %Position of the target with respect to the first primary
            ur1 = Ur1/norm(Ur1);            %Unit vector of the relative position of the target with respect to the primary
            Ur2 = r_t-R(:,2);               %Position of the target with respect to the first primary
            ur2 = Ur2/norm(Ur2);            %Unit vector of the relative position of the target with respect to the primary
            %Evaluate the linear model 
            Sigma = -((mup(1)/norm(Ur1)^3)+(mup(2)/norm(Ur2))^3)*eye(3)+3*((mup(1)/norm(Ur1)^3)*(ur1*ur1.')+(mup(2)/norm(Ur2)^3)*(ur2*ur2.'));
        otherwise 
            error('No valid linear model was selected'); 
    end

    %Linear state model
    A = [zeros(3) eye(3) zeros(3); zeros(3) zeros(3) eye(3); zeros(3) Sigma Omega];  

    %Compute the feedback control law
    [K,~,~] = lqr(A,B,Q,M);
    
    %Compute the control law 
    u = -K*[int; s_r];

    %Compute the integration of the relative motion equations 
    drho = Encke_method(mu, s_t, s_r);                          %Natural vector field flow
    dint = s_r(1:3);                                            %Integrator field flow
    drho = drho + [0; 0; 0; u];                                 %Add control vector
    drho = [drho; dint];
end

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
function [drho] = so_model(mu, s_t, s_r)  
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
function [drho] = th_model(mu, s_t, s_r)  
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