%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 21/05/21
% File: CTR_guidance.m 
% Issue: 0 
% Validated: 20/05/21

%% Chebyshev Trajectory Regression Guidance %%
% This script contains the function to compute the control law by means of the CTRG guidance core.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - string cost_function, for both position, velocity and complete
%           rendezvous: 'Position', 'Velocity', 'State
%         - scalar Tmin, minimum available thrust
%         - scalar Tmax, maximum available thrust
%         - scalar TOF, the time of flight for the rendezvous condition
%         - vector s0, initial conditions of both the target and the
%           relative particle
%         - string core, selecting the solver (linear or nonlinear) to be
%           used
%         - string method, selecting the nonlinear solver to use

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [Cp, Cv, Cg] = CTR_guidance()
end