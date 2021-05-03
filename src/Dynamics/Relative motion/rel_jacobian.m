%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 24/03/21 % 

%% Jacobian %% 
% This script provides a function to compute the jacobian of the relative motion dynamics 
% at a given state. 

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

function [J] = rel_jacobian(mu, s)
    %System parameters 
    mup(1) = 1-mu;                      %Reduced gravitational parameter of the first primary
    mup(2) = mu;                        %Reduced gravitational parameter of the second primary
    R(1:3,1) = [-mu; 0; 0];             %Synodic position of the first primary
    R(1:3,2) = [1-mu; 0; 0];            %Synodic position of the second primary
       
    %State variables 
    r_t = s(1:3);                       %Position of the target 
    rho = s(7:9);                       %Relative position
    
    %Relative position to the primaries 
    r(:,1) = r_t-R(:,1);                %Target position to the first primary
    r(:,2) = r_t-R(:,2);                %Target position to the second primary
    rc(:,1) = rho+r(:,1);               %Chaser position to the first primary
    rc(:,2) = rho+r(:,2);               %Chaser position to the second primary
    
    %Variational equations
    O = zeros(3,3);                     %3 by 3 null matrix
    I = eye(3);                         %3 by 3 identity matrix
    Omega = [0 2 0; -2 0 0; 0 0 0];     %Coriolis dyadic
   
    %Hessian of the potential function
    C = [1; 1; 0];                      %Centrifugal force vector
    H = zeros(3,3);                     %Preallocation of the hessian of the potential 
    for i = 1:size(H,1)
        for j = i:size(H,2)
            %Derivative of the Encke acceleration
            gamma = 0;                                                        %Encke acceleration
            for k = 1:length(mup)
                %Encke's variables
                q = -dot(2*r(:,k)+rho,rho)/norm(rc(:,k))^2;                   %Encke variable
                f = q*(3*(1+q)+q^2)/(1+(1+q)^(3/2));                          %Encke coefficient
                
                %Derivative of the Encke coefficient
                df = -(3*sqrt(1+q)/norm(rc(:,k))^2)*(1-dot(2*(r(:,k)+rho),rho)/norm(rc(:,k))^2)*rc(i,k); 
                
                %Encke acceleration
                if (i == j)
                    gamma = gamma + (mup(k)/norm(r(:,k))^3)*((1+f)+df*rc(j,k));             
                else
                    gamma = gamma + (mup(k)/norm(r(:,k))^3)*(df*rc(j,k));
                end
            end
            
            %Hessian of the potential function
            if (i == j)
                H(i,j) = C(i)-gamma;    %Diagonal terms
            else
                H(i,j) = -gamma;        %Non-diagonal terms
            end
        end
    end
    
    %Symmetry constraint
    for i = 1:size(H,1)
        for j = 1:size(H,2)
            H(j,i) = H(i,j);
        end
    end
    
    %Jacobian of the system 
    J = [O I; H Omega];                 
end