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

function [drho] = opt_model(mu, dt, Sn, t, x, u)
    %Constants of the system 
    mup(1) = 1-mu;                      %Reduced gravitational parameter of the first primary 
    mup(2) = mu;                        %Reduced gravitational parameter of the second primary 
    
    %State variables 
    index = fix(1+t/dt);                %Temporal index
    r_t = Sn(index,1:3).';              %Synodic position of the target
    r_r = x(1:3,:);                     %Synodic relative position 
    v_r = x(4:6,:)+u;                   %Synodic relative velocity 
    
    %Synodic position of the primaries 
    R(1:3,1) = [-mu; 0; 0];             %Synodic position of the first primary
    R(1:3,2) = [1-mu; 0; 0];            %Synodic position of the second primary
    
    %Encke acceleration method
    gamma = [2*v_r(2,:)+r_r(1,:); -2*v_r(1,:)+r_r(2,:); zeros(1,size(r_r,2))];
    for i = 1:length(mup)
        q = -dot(2*(r_t(1:3,:)-R(:,i))+r_r(1:3,:),r_r(1:3,:))/norm(r_t(1:3,:)+r_r(1:3,:)-R(:,i))^2;
        f = q(1,:).*(3*(1+q(1,:))+q(1,:).^2)./(1+(1+q(1,:)).^(3/2));
        gamma = gamma - (mup(i)./norm(r_t(1:3,:)-R(:,i)).^3).*(f(1,:).*(r_t(1:3,:)-R(:,i))+(1+f(1,:)).*r_r(1:3,:));
    end
    
    %Equations of motion 
    drho = [v_r; 
            gamma];
end