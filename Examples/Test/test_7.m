%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/07/21
% File: test_7.m 
% Issue: 0 
% Validated: 

%% Test 7 %%
% This scripts provides a test interface for the rest of the library
% functions. 

% Test 7 is concerned with differential correction algorithms for computation of dynamical 2D tori.

%% Test values and constants
%Set graphical environment 
set_graphics(); 

%Initial conditions
mu = 0.0121505856;                                          %Reduced gravitational parameter of the system (Earth-Moon)
L = libration_points(mu);                                   %System libration points
Az = 100e6;                                                  %Orbit amplitude out of the synodic plane. Play with it!
Ax = 50e6;                                                  %Orbit amplitude in the synodic plane. Play with it! 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ax = dimensionalizer(384400e3, 1, 1, Ax, 'Position', 0);    %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param_halo = [1 Az Ln gamma m];                             %Halo orbit parameters (-1 for southern halo)
param_lyap = [Ax Az 0 0 Ln gamma m];                        %Lyapunov orbit parameters

%Correction parameters 
maxIter = 100;     %Maximum allowed iterations in the differential correction schemes
tol = 1e-5;        %Tolerance 

%% Functions
%Compute seeds
[halo_seed, haloT] = object_seed(mu, param_halo, 'Halo');   %Generate a halo orbit seed
lyapunov_seed = object_seed(mu, param_lyap, 'Lyapunov');    %Generate a Lyapunov orbit seed

%Lyapunov orbit
[lyapunov_orbit, state(1)] = differential_correction('Planar', mu, lyapunov_seed, maxIter, tol);

%Halo orbit 
[halo_orbit, state(2)] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Compute the halo tori 
J = jacobi_constant(mu, halo_orbit.Trajectory(1,1:6).');
[xf, state] = differential_torus('Single shooting energy', mu, halo_orbit, maxIter, tol, 5, J);