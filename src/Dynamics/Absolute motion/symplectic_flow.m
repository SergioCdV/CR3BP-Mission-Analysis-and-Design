%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 15/08/21
% File: symplectic_flow.m 
% Issue: 0 
% Validated: 

%% Symplectic Flow %%
% This function contains the linearized symplectic form of the equations of motion near 
% libration points.

% Inputs: - scalar mu, the sysmtem's gravitational parameter.
%         - scalar L, the libration point near which the motion happens. 
%         - scalar gamma, the characteristic distance associated with L.
%         - scalar t, an integration time span. 
%         - vector x, containing the three-dimensional synodic coordinates
%           to be propagated.

% Output: - S, the propagated synodic coordinates.
%         - z, the vector of propagated Hamiltonian coordinates.

% Methods: linearized second order Hamiltonian flow in symplectic
%          coordinates

% New version: variational equations 

function [S, z] = symplectic_flow(mu, L, gamma, t, x)
    %Compute the eigenvalues of the problem 
    m = 6;                                                  %Dimension of the system 
    c = legendre_coefficients(mu, L, gamma, 2);             %Second order Legendre coefficient
    c2 = c(end);                                            %Second order Legendre coefficient
    alpha(1) = (c2-2+sqrt(9*c2^2-8*c2))/2;                  %Characteristic root
    alpha(2) = (c2-2-sqrt(9*c2^2-8*c2))/2;                  %Characteristic root
    
    lambda = sqrt(alpha(1));                                %Hyperbolic eigenvalue
    omegap = sqrt(-alpha(2));                               %Center manifold eigenvalue 
    omegav = sqrt(c2);                                      %Center manifold eigenvalue
    
    %Compute the linearizing matrix
    s(1) = 2*lambda*((4+3*c2)*lambda^2+4+5*c2-6*c2^2);      %Scaling factor
    s(2) = omegap*((4+3*c2)*omegap^2-4-5*c2+6*c2^2);        %Scaling factor
    s = sqrt(s);
    C = zeros(m,m);
    C(1,1) = 2*lambda/s(1); 
    C(2,1) = (lambda^2-2*c2-1)/s(1);
    C(4,1) = (lambda^2+2*c2+1)/s(1);
    C(5,1) = (lambda^3+(1-2*c2)*lambda)/s(1);
    C(2,2) = (-omegap^2-2*c2-1)/s(2);
    C(4,2) = (-omegap^2+2*c2+1)/s(2);
    C(3,3) = 1/sqrt(omegav);   
    C(1,4) = -2*lambda/s(1); 
    C(2,4) = (lambda^2-2*c2-1)/s(1);
    C(4,4) = (lambda^2+2*c2+1)/s(1);
    C(5,4) = (-lambda^3-(1-2*c2)*lambda)/s(1);
    C(1,5) = 2*omegap/s(2);
    C(5,5) = (-omegap^3+(1-2*c2)*omegap)/s(2);
    C(6,6) = sqrt(omegav);
   
    %Preallocation 
    z = zeros(m,length(t));     %Hamiltonian state vector 
    S = zeros(length(t),m);     %Synodic state vector
    
    %Hamiltonian and symplectic state vectors 
    z0(1,1) = x(1);             %Synodic x coordinate
    z0(2,1) = x(2);             %Synodic y coordinate
    z0(3,1) = x(3);             %Synodic z coordinate
    z0(4,1) = x(4)-x(2);        %Synodic x momenta
    z0(5,1) = x(5)+x(1);        %Synodic y momenta
    z0(6,1) = x(6);             %Synodic z momenta
    z0 = C^(-1)*z0;             %Generalized symplectic state vector 
    
    %Propagation of the symplectic coordinates
    z(1,:) = z0(1)*exp(lambda*t);
    z(2,:) = real((z0(2)+1i*z0(5))*exp(-1i*omegap*t));
    z(3,:) = real((z0(3)+1i*z0(6))*exp(-1i*omegav*t));
    z(4,:) = z0(4)*exp(-lambda*t);
    z(5,:) = imag((z0(2)+1i*z0(5))*exp(-1i*omegap*t));
    z(6,:) = imag((z0(3)+1i*z0(6))*exp(-1i*omegav*t));
    
    for i = 1:length(t)
        z(:,i) = C*z(:,i);
    end

    %Synodic coordinates re-mapping
    S(:,1:3) = z(1:3,:).';          %Synodic position coordinates
    S(:,4) = (z(4,:)+z(2,:)).';     %Synodic x velocity
    S(:,5) = (z(5,:)-z(1,:)).';     %Synodic y velocity
    S(:,6) = z(6,:).';              %Synodic z velocity
end