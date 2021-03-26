%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/03/21
% File: cr3bp_jacobian.m 
% Issue: 0 
% Validated: 

%% CR3BP Jacobian %%
% This function contains the jacobian of the vector field of the CR3BP. It accounts for a infinitesimal mass
% moving in the normalized, non dimensional synodic frame define by the two primaries, which
% are assumed to be in the same plane and in circular orbits.

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - vector s, containing in an Mx1 array the phase space vector, with M 
%           the phase space dimension.

% Outputs: - matrix Jacobian, containing the Jacobian evaluated at the phase space vector.

% Methods: non-dimensional CR3BP dynamics in the synodic frame. 

% New versions: 

function [Jacobian] = cr3bp_jacobian(mu, s) 
    %Define the phase space vector
    x = s(1);               %Synodyc x coordinate
    y = s(2);               %Synodyc y coordinate 
    z = s(3);               %Synodyc z coordinate 
    
    %Relevant system parameters
    mu1 = 1-mu;             %First primary normalized position
    mu2 = mu;               %Second primary normalized position
    r1 = [(x+mu2); y; z];   %Relative position vector to the first primary
    r2 = [(x-mu1); y; z];   %Relative position vector to the secondary primary
    R1 = norm(r1);          %Distance to the first primary
    R2 = norm(r2);          %Distance to the secondary primary
    
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
    Jacobian = [O I; G K];                             %Jacobian of the system   
end