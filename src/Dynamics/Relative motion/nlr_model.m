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

function [ds] = nlr_model(mu, direction, flagVar, method_ID, t, s)
    %State variables 
    s_t = s(1:6);       %State of the target
    s_r = s(7:12);      %State of the relative particle
    
    %Equations of motion of the target
    ds_t = cr3bp_equations(mu, direction, flagVar, t, s_t);       %Target equations of motion
    
    %Equations of motion of the relative state 
    switch (method_ID)
        case 'Encke'
            drho = Encke_method(mu, s_t, s_r);                    %Relative motion equations
        case 'Encke V'
            drho = EnckeV_method(mu, s);                          %Relative motion equations
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
function [drho] = EnckeV_method(mu, s)
    %System parameters 
    m = 6;                          %Phase space dimension
    mup(1) = 1-mu;                  %Reduced gravitational parameter of the first primary
    mup(2) = mu;                    %Reduced gravitational parameter of the second primary
    R(1:3,1) = [-mu; 0; 0];         %Synodic position of the first primary
    R(1:3,2) = [1-mu; 0; 0];        %Synodic position of the second primary
       
    %State variables 
    r_t = s(1:3);       %Position of the target 
    s_t = s(1:6);       %State of the target
    rho = s(7:9);       %Relative position
    s_r = s(7:12);      %Relative state
    ds = s(13:end);     %Phase space error vector 
    
    %Relative position to the primaries 
    r(:,1) = r_t-R(:,1);    %Target position to the first primary
    r(:,2) = r_t-R(:,2);    %Target position to the second primary
    rc(:,1) = rho+r(:,1);   %Chaser position to the first primary
    rc(:,2) = rho+r(:,2);   %Chaser position to the second primary

    %Compute the integration of the relative motion equations 
    drho = Encke_method(mu, s_t, s_r);
    
    %Variational equations
    ds = reshape(ds, [m m]);            %Reshape the covariance vector in the matrix form 
    O = zeros(3,3);                     %3 by 3 null matrix
    I = eye(3);                         %3 by 3 identity matrix
    Omega = [0 2 0; -2 0 0; 0 0 0];     %Coriolis dyadic
   
    %Hessian of the potential function
    C = [1; 1; 0];                      %Centrifugal force vector
    H = zeros(3,3);                     %Preallocation of the hessian of the potential 
    for i = 1:size(H,1)
        for j = i:size(H,2)
            %Derivative of the Encke acceleration
            gamma = 0;                  %Encke acceleration
            for k = 1:length(mup)
                q = -dot(2*r(:,k)+rho,rho)/norm(rc(:,k))^2;                   %Encke variable
                f = q*(3*(1+q)+q^2)/(1+(1+q)^(3/2));                          %Encke coefficient
                %Derivative of the Encke coefficient
                df = -(3*sqrt(1+q)/norm(rc(:,k))^2)*(1-dot(r(:,k),rho)/norm(rc(:,k))^2)*rc(j,k);                             
                %Encke acceleration
                if (i == j)
                    gamma = gamma + (mup(k)/norm(r(:,k))^3)*((1+f)+df*rc(j,k));             
                else
                    gamma = gamma + (mup(k)/norm(r(:,k))^3)*(df*rc(j,k));
                end
            end
            
            %Hessian of the potential function
            if (i == j)
                H(i,j) = C(i)-gamma;
            else
                H(i,j) = -gamma;
            end
            
            H(i,j) = H(j,i);            %Symmetry constraint
        end
    end
    
    Jac = [O I; H Omega];               %Jacobian of the system 
    ds = Jac*ds;                        %Variational equations 
    ds = reshape(ds, [m^2 1]);          %Variational vector field 
    
    %Final vector field 
    drho = [drho; ds]; 
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
    Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];                                    %Linear term
    Sigma3 = [0; 0; 0; 3*c3*x^2-(3/2)*c3*(y^2+z^2); -3*c3*x*y; -3*c3*x*z];      %Second order term
    
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
    Sigma3 = [0; 0; 0; 3*c3*x^2-(3/2)*c3*(y^2+z^2); -3*c3*x*y; -3*c3*x*z];      %Second order term
    Sigma4 = [0; 0; 0; ...
              (5/2)*c4*x^3-(30/4)*c4*(y^2+z^2)*x+(3/2)*c4*(x^2+y^2+z^2)*x; ...
             -(30/4)*c4*x^2*y+(3/2)*c4*(x^2+y^2+z^2)*y; ...
             -(30/4)*c4*x^2*z+(3/2)*c4*(x^2+y^2+z^2)*z];                        %Third order term
    
    %Equations of motion 
    drho = [O I; Sigma -2*Omega]*s_r + Sigma3 + Sigma4;
end