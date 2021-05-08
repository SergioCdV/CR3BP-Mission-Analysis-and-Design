%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/05/21
% File: LQR_control.m 
% Issue: 0 
% Validated: 08/05/21

%% State Dependt Ricatti Equation Control %%
% This script contains the function to compute the control law by means of an SDRE controller.

% Inputs: - string model, selecting the linear model to compute the linear state
%           transition matrix of the system
%         - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - array Sn, the target orbit evolution 
%         - scalar Ln, the libration point number. It may be left 0 if the
%           target orbit is not librating around any Lagrange point
%         - scalar gamma, the relative distance of the Ln point to the
%           nearest primary. Again, it may be left as 0 if needed

% Output: - vector u, the computed control law

% New versions: 
