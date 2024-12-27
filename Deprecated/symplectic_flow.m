%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 15/08/21
% File: symplectic_flow.m 
% Issue: 0 
% Validated: 

%% Symplectic Flow %%
% This function contains the linearized symplectic form of the equations of motion near 
% libration points

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
    % Compute the eigenvalues of the problem 
    c = legendre_coefficients(mu, L, gamma, 2);             % Second order Legendre coefficient
    c2 = c(end);                                            % Second order Legendre coefficient
    alpha(1) = (c2-2+sqrt(9*c2^2-8*c2))/2;                  % Characteristic root
    alpha(2) = (c2-2-sqrt(9*c2^2-8*c2))/2;                  % Characteristic root
    
    lambda = sqrt(alpha(1));                                % Hyperbolic eigenvalue
    omegap = sqrt(-alpha(2));                               % Center manifold eigenvalue 
    omegav = sqrt(c2);                                      % Center manifold eigenvalue
    
    % Reference frame transformation
    [z0, C] = symplectic_frame(mu, L, gamma, x);

    % Propagation of the symplectic coordinates
    z(1,:) = z0(1)*exp(lambda*t);
    z(2,:) = real((z0(2)+1i*z0(5))*exp(-1i*omegap*t));
    z(3,:) = real((z0(3)+1i*z0(6))*exp(-1i*omegav*t));
    z(4,:) = z0(4)*exp(-lambda*t);
    z(5,:) = imag((z0(2)+1i*z0(5))*exp(-1i*omegap*t));
    z(6,:) = imag((z0(3)+1i*z0(6))*exp(-1i*omegav*t));
   
    % Synodic coordinates re-mapping
    z = C*z;
    S(:,1:3) = z(1:3,:).';          % Synodic position coordinates
    S(:,4) = (z(4,:)+z(2,:)).';     % Synodic x velocity
    S(:,5) = (z(5,:)-z(1,:)).';     % Synodic y velocity
    S(:,6) = z(6,:).';              % Synodic z velocity
end