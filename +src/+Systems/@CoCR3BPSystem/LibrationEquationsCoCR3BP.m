%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: LibrationEquationsCoCR3BP.m 
% Issue: 0 
% Validated: 

%% Linear Equations of the Co-orbital CR3BP Dynamics around a libration point %%
% This function contains the description of the co-orbital CR3BP dynamics vector field

% Inputs: 
% Outputs: - vector ds, the differential vector field

% New versions: 

function [ds] = LibrationEquationsCoCR3BP(t, j, s, u, params)
    % State variables 
    r_t = s(1:3);                               % Position vector of the target
    s_r = s(7:12,:);                            % Relative synodic state vector
        
    % Relative Legendre coefficient c2           
    mu = params(1);                                                      % Gravitational parameter of the system
    cn = src.Systems.CoCR3BPSystem.CoLegendreCoefficients(mu, r_t, 2);   % Relative Legendre coefficients 
    c2 = cn(2);                                                          % First order relative Legendre coefficient
    
    % Relative acceleration (non inertial)
    O = zeros(3);                                       % 3 by 3 null matrix
    I = eye(3);                                         % 3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];                     % Hat map dyadic of the angular velocity for the synodice reference frame
    Sigma = [1 + 2 * c2 0 0; 0 1 - c2 0; 0 0 -c2];      % Gravity acceleration
    
    % Equations of motion 
    ds = [O I; Sigma 2*Omega] * s_r;
 
    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;
end