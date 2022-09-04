%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/11/21
% File: synmplectic_frame.m 
% Issue: 0 
% Validated: 

%% Symplectic frame %%
% For a given gravitational parameter mu and a space state vector and a given libration point, for the CR3BP dynamics, this function
% computes the associated symplectic state vector

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - scalar L, the identifier of the orbit libration point
%         - scalar gamma, the distance of the orbit libration point to one
%           of the primaries
%         - vector s, phase space vector

% Outputs: - the vector S, the state vector in the symplectic coordinate
%            frame
%          - matrix C, the symplectic-synodic canonical transformation
%            matrix

% New versions: 

function [S, C] = symplectic_frame(mu, L, gamma, s)
    % Compute the eigenvalues of the problem 
    m = 6;                                                  % Dimension of the system 
    c = legendre_coefficients(mu, L, gamma, 2);             % Second order Legendre coefficient
    c2 = c(end);                                            % Second order Legendre coefficient
    alpha(1) = (c2-2+sqrt(9*c2^2-8*c2))/2;                  % Characteristic root
    alpha(2) = (c2-2-sqrt(9*c2^2-8*c2))/2;                  % Characteristic root
    
    lambda = sqrt(alpha(1));                                % Hyperbolic unstable eigenvalue
    omegap = sqrt(-alpha(2));                               % Hyperbolic stable eigenvalue
    omegav = sqrt(c2);                                      % Center manifold eigenvalue
    
    % Compute the canonical transformation matrix
    s(1) = 2*lambda*((4+3*c2)*lambda^2+4+5*c2-6*c2^2);      % Scaling factor
    s(2) = omegap*((4+3*c2)*omegap^2-4-5*c2+6*c2^2);        % Scaling factor
    s = sqrt(s);                                            % Scaling factor

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
    
    % Hamiltonian and symplectic state vectors 
    z(1,1) = s(1);             % Synodic x coordinate
    z(2,1) = s(2);             % Synodic y coordinate
    z(3,1) = s(3);             % Synodic z coordinate
    z(4,1) = s(4)-s(2);        % Synodic x momenta
    z(5,1) = s(5)+s(1);        % Synodic y momenta
    z(6,1) = s(6);             % Synodic z momenta
    S = C^(-1)*z;              % Generalized symplectic state vector   
end