%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Initial approximation %%
% Function to estimate the initial time of flight, control points and curve approximation

% Inputs: - class Problem, defining the problem of interest
%         - string basis, specifying the polynomial collacation basis
%         - vector tau, the collocation points to be used 

% Outputs: - vector betaapp, the initial estimation of the optimization
%            extra variables
%          - scalar tfapp, the initial initial time 
%          - scalar tfapp, the initial initial time of flight
%          - array Papp, the initial estimation of the boundary control
%            points
%          - array Capp, the initial estimation of the spacecraft state vector

function [betaapp, t0app, tfapp, Papp, Capp] = initial_approximation(Problem, basis, tau)
    % Constants 
    L = Problem.DerDeg;         % Order of the dynamics (maximum derivative order)

    % Initial guess 
    [betaapp, t0app, tfapp] = Problem.InitialGuess(Problem.Params, Problem.initial, Problem.final);

    % Generate the polynomial basis
    n_init = repmat(3, [1 Problem.StateDim]).';
    Bapp = state_basis(n_init, L, basis, tau);

    % Initial estimate of control points (using the non-orthonormal boundary conditions)
    Papp = zeros(Problem.StateDim, max(n_init)+1);  
    Papp = boundary_conditions(Problem, betaapp, t0app, tfapp, Bapp, basis, n_init, Papp);

    % State vector approximation as a function of time
    Capp = evaluate_state(Papp, Bapp, n_init, L);
end