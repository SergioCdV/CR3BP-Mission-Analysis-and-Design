%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
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
%         - scalar t, a reference epoch. 
%         - vector s, containing in an Nx1 array the phase space vector,
%           possibly augmented with the state transition matrix at time t. 

% Outputs: - vector dr, the differential vector field, which will include
%            the phase space trajectory.

% Methods: . 

% New versions: include the first variations of the vector field.

function [ds] = fulrel_motion(mu, direction, flagVar, Encke_flag, t, s)
    %State variables 
    s_t = s(1:6);       %State of the target
    rho = s(7:12);      %State of the chaser
    
    %Equations of motion of the target
    ds_t = cr3bp_equations(mu, direction, flagVar, t, s_t);        %Target equations of motion
    
    %Equations of motion of the relative state 
    if (Encke_flag)
        drho = Encke_method(mu, s_t, rho);                         %Relative motion equations
    else
        drho = relative_motion(mu, s_t, rho);                      %Relative motion equations
    end
    
    %Vector field 
    ds = [ds_t; drho];
end

%% Auxiliary functions 
%Relative motion equations
function [drho] = relative_motion(mu, s_t, s_r)
    %Constants of the system 
    mu1 = 1-mu;             %Reduced gravitational parameter of the first primary 
    mu2 = mu;               %Reduced gravitational parameter of the second primary 
    
    %State variables 
    r_t = s_t(1:3);         %Synodic position of the target
    r_r = s_r(1:3);         %Synodic relative position 
    v_r = s_r(4:6);         %Synodic relative velotice 
    
    %Synodic position of the primaries 
    R1 = [-mu; 0; 0];       %Synodic position of the first primary
    R2 = [1-mu; 0; 0];      %Synodic position of the second primary
    
    %Equations of motion 
    drho = [v_r; 
            r_r(1)+2*v_r(2)-mu1*((r_t(1)+mu)/norm(r_t-R1)^3-(r_t(1)+r_r(1)+mu)/norm(r_t-R1+r_r)^3)-mu2*((r_t(1)-(1-mu))/norm(r_t-R2)^3-(r_t(1)+r_r(1)-(1-mu))/norm(r_t-R2+r_r)^3); 
            r_r(2)-2*v_r(1)-mu1*(r_t(2)/norm(r_t-R1)^3-(r_t(2)+r_r(2))/norm(r_t-R1+r_r)^3)-mu2*(r_t(2)/norm(r_t-R2)^3-(r_t(2)+r_r(2))/norm(r_t-R2+r_r)^3); 
            -mu1*(r_t(3)/norm(r_t-R1)^3-(r_t(3)+r_r(3))/norm(r_t-R1+r_r)^3)-mu2*(r_t(3)/norm(r_t-R2)^3-(r_t(3)+r_r(3))/norm(r_t-R2+r_r)^3)];
end

function [drho] = Encke_method(mu, s_t, s_r)
    %Constants of the system 
    mu_r(1) = 1-mu;                 %Reduced gravitational parameter of the first primary 
    mu_r(2) = mu;                   %Reduced gravitational parameter of the second primary 
    
    %State variables 
    r_t = s_t(1:3);               %Synodic position of the target
    r_r = s_r(1:3);               %Synodic relative position 
    v_r = s_r(4:6);               %Synodic relative velotice 
    
    %Synodic position of the primaries 
    R(1:3,1) = [-mu; 0; 0];       %Synodic position of the first primary
    R(1:3,2) = [1-mu; 0; 0];      %Synodic position of the second primary
    
    %Encke acceleration method
    gamma = [2*v_r(2)+r_r(1); -2*v_r(1)+r_r(2); 0];
    for i = 1:length(mu_r)
        q = -dot(2*(r_t-R(:,i))+r_r,r_r)/norm(r_t+r_r-R(:,i))^2;
        f = q*(3*(1+q)+q^2)/(1+(1+q)^(3/2));
        gamma = gamma + (mu_r(i)/norm(r_t-R(:,i))^3)*(f*(r_t-R(:,i))+(1+f)*r_r);
    end
    
    %Equations of motion 
    drho = [v_r; gamma];
end