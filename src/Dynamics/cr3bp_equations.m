%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: cr3bp_equations.m 
% Issue: 0 
% Validated: 

%% CR3BP Dynamics %%
% This function contains the 42-DOF translational model used within the rest of
% the CR3BP General Library scripts. It accounts for a infinitesimal mass
% moving in the normalized, non dimensional synodic frame define by the two primaries, which
% are assumed to be in the same plane and in circular orbits. It also
% contains the integration of the first variational equations of the flow.

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
%            the phase space trajectory and the STM integrated state when
%            flagVar is true.

% Methods: non-dimensional CR3BP dynamics in the synodic frame. 

% New versions: 

function [dr] = cr3bp_equations(mu, direction, flagVar, t, s)
    %Constants 
    n = 6;                  %Phase space dimension 
    
    %Define the initial phase space vector
    x = s(1);               %Synodyc x coordinate
    y = s(2);               %Synodyc y coordinate 
    z = s(3);               %Synodyc z coordinate 
    V = s(4:6);             %Synodic velocity vector
    
    %Relevant system parameters
    mu1 = 1-mu;             %First primary normalized position
    mu2 = mu;               %Second primary normalized position
    r1 = [(x+mu2); y; z];   %Relative position vector to the first primary
    r2 = [(x-mu1); y; z];   %Relative position vector to the secondary primary
    R1 = norm(r1);          %Distance to the first primary
    R2 = norm(r2);          %Distance to the secondary primary
    
    %Compute the time flow of the system (depending on the time direction)
    if (direction == -1)
        inAcc = [x-2*V(2); y+2*V(1); 0];                    %Inertial acceleration
    else
        inAcc = [x+2*V(2); y-2*V(1); 0];                    %Inertial acceleration
    end
    F = [V; inAcc-mu1/R1^3*r1-mu2/R2^3*r2];                 %Time flow of the system
    
    %Compute the variational equations if needed
    if (flagVar)
        %Compute the initial STM
        phi = reshape(s(n+1:end), [n n]);
        
        %First variations of the augmented potential function (Hessian of the potential)
        Ux = 1-(mu1/R1^3)*(1-3*((x+mu2)/R1)^2)-(mu2/R2^3)*(1-3*((x-mu1)/R2)^2); 
        Uy = 1-(mu1/R1^3)*(1-3*(y/R1)^2)-(mu2/R2^3)*(1-3*(y/R2)^2);
        Uz = -(mu1/R1^3)*(1-3*(z/R1)^2)-(mu2/R2^3)*(1-3*(z/R2)^2);
        Uxy = 3*y*((mu1/R1^5)*(x+mu2)+(mu2/R2^5)*(x-mu1));
        Uyx = Uxy;
        Uxz = 3*z*((mu1/R1^5)*(x+mu2)+(mu2/R2^5)*(x-mu1));
        Uzx = Uxz; 
        Uyz = 3*y*((mu1/R1^5)*z+(mu2/R2^5)*z);
        Uzy = Uyz;
        
        %Compute the first variational equations evaluated at the reference
        O = zeros(3,3); 
        I = eye(3);
        G = [Ux Uxy Uxz; Uyx Uy Uyz; Uzx Uzy Uz];
        K = [0 2 0; -2 0 0; 0 0 0];
        Jacob = [O I; G K];                             %Jacobian of the system   
        dphi = Jacob*phi;                               %Variational equations
        dphi = reshape(dphi, [n^2 1]); 
        
        %Update the differential configuration space vector
        dr = [F; dphi];
        
        %Reverse the flow for backward integration
        if (direction == -1)
            dr = -dr;
        end
        
    else
        %Update the differential configuration space vector
        dr = F;  
        
        %Reverse the flow for backward integration
        if (direction == -1)
            dr = -dr;
        end
    end
end