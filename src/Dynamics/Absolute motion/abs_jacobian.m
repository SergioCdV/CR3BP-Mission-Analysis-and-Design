%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/03/21
% File: abs_jacobian.m 
% Issue: 0 
% Validated: 

%% Absolute Jacobian %%
% This function contains the jacobian of the vector field of the CR3BP. It accounts for a infinitesimal mass
% moving in the normalized, non dimensional synodic frame define by the two primaries, which
% are assumed to be in the same plane and in circular orbits.

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - vector s, containing in an Mx1 array the phase space vector, with M 
%           the phase space dimension.

% Outputs: - matrix Jacobian, containing the Jacobian evaluated at the phase space vector.

% Methods: non-dimensional CR3BP dynamics in the synodic frame. 

% New versions: 

function [Jacobian] = abs_jacobian(mu, s) 
    %Define the phase space vector
    x = s(1);                       %Synodyc x coordinate
    y = s(2);                       %Synodyc y coordinate 
    z = s(3);                       %Synodyc z coordinate 
    
    %Relevant system parameters
    mup(1) = 1-mu;                  %First primary normalized position
    mup(2) = mu;                    %Second primary normalized position
    r(:,1) = [x+mup(2); y; z];      %Relative position vector to the first primary
    r(:,2) = [x-mup(1); y; z];      %Relative position vector to the secondary primary
    R(1) = norm(r(:,1));            %Distance to the first primary
    R(2) = norm(r(:,2));            %Distance to the secondary primary
    
    %First variations of the augmented potential function (Hessian of the potential)
    G(1,1) = 1-(mup(1)/R(1)^3)*(1-3*((x+mup(2))/R(1))^2)-(mup(2)/R(2)^3)*(1-3*((x-mup(1))/R(2))^2); 
    G(2,2) = 1-(mup(1)/R(1)^3)*(1-3*(y/R(1))^2)-(mup(2)/R(2)^3)*(1-3*(y/R(2))^2);
    G(3,3) = -(mup(1)/R(1)^3)*(1-3*(z/R(1))^2)-(mup(2)/R(2)^3)*(1-3*(z/R(2))^2);
    G(1,2) = 3*y*((mup(1)/R(1)^5)*(x+mup(2))+(mup(2)/R(2)^5)*(x-mup(1)));
    G(2,1) = G(1,2);
    G(1,3) = 3*z*((mup(1)/R(1)^5)*(x+mup(2))+(mup(2)/R(2)^5)*(x-mup(1)));
    G(3,1) = G(1,3); 
    G(2,3) = 3*y*((mup(1)/R(1)^5)*z+(mup(2)/R(2)^5)*z);
    G(3,2) = G(2,3);

    %Compute the first variational equations evaluated at the reference
    O = zeros(3,3);                    %Null matrix
    I = eye(3);                        %Identity matrix
    K = [0 2 0; -2 0 0; 0 0 0];        %Coriolis Dyadic
    Jacobian = [O I; G K];             %Jacobian of the system   
end