%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: RichardsonEquationsCoCR3BP.m 
% Issue: 0 
% Validated: 

%% Linear Equations of the Co-orbital CR3BP Dynamics around a libration point %%
% This function contains the description of the co-orbital CR3BP dynamics vector field

% Inputs: 
% Outputs: - vector ds, the differential vector field

% New versions: 

function [ds] = RichardsonEquationsCoCR3BP(t, j, s, u, params)
    % Define the initial phase space vector
    s_r = s(7:12,:);                                     % Relative synodic state vector
    
    % Constants of the system 
    c2 = params(end);

    % Relative acceleration vector field
    O = zeros(3);                                       % 3 by 3 null matrix
    I = eye(3);                                         % 3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];                     % Hat map dyadic of the angular velocity for the synodice reference frame
    Sigma = [1 + 2 * c2 0 0; 0 1 - c2 0; 0 0 -c2];      % Gravity acceleration
    A = [O I; Sigma 2 * Omega];                         % Constant state matrix 

    % Relative acceleration (linear order term)
    ds = A * s_r;
 
    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;
end