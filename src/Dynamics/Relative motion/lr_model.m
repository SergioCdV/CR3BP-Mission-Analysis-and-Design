%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/20
% File: rel_motion.m 
% Issue: 0 
% Validated: 

%% Linear relative motion model in the CR3BP %%
% This function contains several linear relative motion models between two punctual particles in the 
% circular restricted three body problem. It accounts for the two masses moving in the normalized, 
% non dimensional synodic frame define by the two primaries, which are assumed to be in the same 
% plane and in circular orbits. 

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - scalar direction (in binary format, 1 or -1), indicating the
%           time integration direction: 1 for forward integration, -1 for
%           backward integration.
%         - boolean flagVar, true for dyanmics and STM integration, 
%           false for only dynamical integration.
%         - scalar t, a reference epoch. 
%         - vector s, containing in an Nx1 array the phase space vector,
%           possibly augmented with the state transition matrix at time t. 

% Outputs: - vector dr, the differential vector field, which will include
%            the phase space trajectory.

% Methods: different linear models of relative motion. 

% New versions: include the first variations of the vector field.

function [ds] = lr_model(mu, cn, direction, flagVar, model, t, s)
    %State variables 
    s_t = s(1:6);       %State of the target
    s_r = s(7:12);      %State of the chaser
    
    %Equations of motion of the target
    ds_t = cr3bp_equations(mu, direction, flagVar, t, s_t);        %Target equations of motion
    
    %Equations of motion of the relative state 
    switch (model)
        case 'Target'
            drho = target_centered(mu, s_t, s_r);                  %Relative motion equations
        case 'Fixed libration'
            drho = librationf_centered(cn, s_r);                   %Relative motion equations
        case 'Moving libration'
            drho = librationm_centered(mu, s_t, s_r);              %Relative motion equations
        otherwise
            drho = [];
            disp('No valid model was chosen');
    end
    
    %Vector field 
    ds = [ds_t; drho];
end

%% Auxiliary functions 
%Relative motion equations linearized with respect to the target
function [drho] = target_centered(mu, s_t, s_r)
    %Constants of the system 
    mup(1) = 1-mu;             %Reduced gravitational parameter of the first primary 
    mup(2) = mu;               %Reduced gravitational parameter of the second primary 
    
    %State variables 
    r_t = s_t(1:3);            %Synodic position of the target
    
    %Synodic position of the primaries 
    R(:,1) = [-mu; 0; 0];      %Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];     %Synodic position of the second primary
    
    %Relative position between the primaries and the target 
    Ur1 = r_t-R(:,1);                       ¡%Position of the target with respect to the first primary
    ur1 = Ur1/norm(Ur1);                    %Unit vector of the relative position of the target with respect to the first primary
    Ur2 = r_t-R(:,2);                       %Position of the target with respect to the first primary
    ur2 = Ur2/norm(Ur2);                    %Unit vector of the relative position of the target with respect to the second primary
    
    %Relative acceleration (non inertial)
    O = zeros(3,3);                         %3 by 3 null matrix
    I = eye(3);                             %3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];         %Hat map dyadic of the angular velocity for the synodice reference frame
    
    %Gravity acceleration
    Sigma = -((mup(1)/norm(Ur1)^3)+(mup(2)/norm(Ur2))^3)*eye(3)+3*((mup(1)/norm(Ur1)^3)*(ur1*ur1.')+(mup(2)/norm(Ur2)^3)*(ur2*ur2.'));
    
    %State matrix 
    A = [O I; Sigma-Omega*Omega -2*Omega];
    
    %Equations of motion 
    drho = A*s_r;
end

%Relative motion equations linearized with respect to the global libration point
function [drho] = librationf_centered(cn, s_r)    
    %Legendre coefficient c2           
    c2 = cn(2);                             %Second Legendre coefficient
    
    %Relative acceleration (non inertial)
    O = zeros(3,3);                         %3 by 3 null matrix
    I = eye(3);                             %3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];         %Hat map dyadic of the angular velocity for the synodice reference frame
    
    %Gravity acceleration
    Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];
    
    %Equations of motion 
    drho = [O I; Sigma -2*Omega]*s_r;
end

%Relative motion equations linearized with respect to the local libration point
function [drho] = librationm_centered(mu, s_t, s_r)  
    %State variables 
    r_t = s_t(1:3);                             %Position vector of the target
        
    %Relative Legendre coefficient c2           
    cn = relegendre_coefficients(mu, r_t, 2);   %Relative Legendre coefficients 
    c2 = cn(2);                                 %First order relative Legendre coefficient
    
    %Relative acceleration (non inertial)
    O = zeros(3,3);                             %3 by 3 null matrix
    I = eye(3);                                 %3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];             %Hat map dyadic of the angular velocity for the synodice reference frame
    
    %Gravity acceleration
    Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];
    
    %Equations of motion 
    drho = [O I; Sigma -2*Omega]*s_r;
end

