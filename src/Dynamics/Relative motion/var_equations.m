%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 03/10/22
% File: var_equations.m 
% Issue: 0 
% Validated: 

%% Linear variational relative motion equations in the CR3BP %%
% This function contains the full nonlinear vectorized relative motion between two punctual particles in the 
% CR3BP. It accounts for the two masses moving in the normalized, non dimensional synodic frame 
% define by the two primaries, which are assumed to be in the same plane and in circular orbits

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - scalar direction (in binary format, 1 or -1), indicating the
%           time integration direction: 1 for forward integration, -1 for
%           backward integration
%         - boolean flagVar, true for dyanmics and STM integration, 
%           false for only dynamical integration
%         - boolean relFlagVar, true for relative dyanmics and STM integration, 
%           false for only dynamical integration
%         - string method_ID, identifying which integration method to use
%         - scalar t, a reference epoch
%         - vector s, containing in an Nx1 array the phase space vector,
%           possibly augmented with the state transition matrix at time t
%         - cell array varargin, to include GNC requirements on the motion
%           of the spacecraft

% Outputs: - vector ds, the differential vector field, which will include
%            the phase space trajectory

% Methods: 

% New versions: vectorize all models and options

% Full nonlinear relative motion equations via Encke's method
function [dphi] = var_equations(mu, t, s)
    % Constants of the system 
    n = 6;                        % Phase space dimension

    s_t = s(1:n,:);               % Target state
    s_r = s(n+1:end,:);           % Chaser relative state

    if (mod(size(s_r,1),n+n/2) == 0)
        n = n+n/2;
    end

    % Variational equations
    dphi = zeros(n^2,size(s_r,2));          % STM preallocation
    J = zeros(n,n) ;                        % Preallocation of the Jacobian matrix
    if (n == 9)
        J(7:9,1:3) = eye(3);
    end

    for i = 1:size(s_r,2)
        Phi = reshape(s_r(n+1:end,i), [n n]);                        % State Transition Matrix
        J(1:6,1:6) = rel_jacobian(mu, [s_t(:,i); s_r(1:3,i)]);       % Relative Jacobian matrix
        aux = J*Phi;                                                 % New State Transition Matrix
        dphi(:,i) = reshape(aux, [n^2 1]);                           % First linear variational vector field
    end
end
