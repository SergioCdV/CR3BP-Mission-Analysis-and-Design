%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 28/11/22 % 

%% Numerical error %% 
% This script provides a script to test different nonlinear models for relative motion in the CR3BP. 

% The relative motion of two spacecraft in a halo orbit around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion % 

% Units are non-dimensional and solutions are expressed in the synodic reference frame as defined by Howell, 1984.

%% Set up %%
close all; 
clear; 
clc; 

% Set up graphics 
set_graphics();

% Integration tolerances
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Constants and initial data %% 
% Time span 
dt = 1e-3;                          % Time step
tmax = pi;                          % Maximum time of integration (corresponding to a synodic period)
tspan = 0:dt:tmax;                  % Integration time span

% CR3BP constants 
mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
L = libration_points(mu);           % System libration points
Lem = 384400e3;                     % Mean distance from the Earth to the Moon

% Differential corrector set up
nodes = 10;                         % Number of nodes for the multiple shooting corrector
maxIter = 20;                       % Maximum number of iterations
tol = 1e-10;                        % Differential corrector tolerance

%% Nominal orbits computation %%
% Halo characteristics 
Az = 20e6;                                                  % Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         % Normalize distances for the E-M system
Ln = 1;                                                     % Orbits around L1
gamma = L(end,Ln);                                          % Li distance to the second primary
m = 1;                                                      % Number of periods to compute

% Compute a halo orbit seed 
halo_param = [1 Az Ln gamma m];                             % Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  % Generate a halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[halo_orbit, state] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

% Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                                     % Bifucartion tolerance on the stability index
num = 2;                                                            % Number of orbits to continuate
method = 'SPC';                                                     % Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                        % Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, period};                              % Object and characteristics to continuate
corrector = 'Plane Symmetric';                                      % Differential corrector method
direction = 1;                                                      % Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                                 % General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(end,:), maxIter, tol);

%% Models comparison %% 
r_t0 = halo_orbit.Trajectory(1,1:6);                        % Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);                      % Initial chaser conditions 
rho0 = r_c0-r_t0;                                           % Initial relative conditions
s0 = [r_t0 rho0];                                           % Initial conditions of the target and the relative state

% Integration of the double-precision relative models
[~, S_c] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, r_c0, options);
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
[t, Sn] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Full nonlinear', t, s), tspan, s0, options);

% Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                 % Reconstructed chaser motion via Encke method
S_rcn = Sn(:,1:6)+Sn(:,7:12);              % Reconstructed chaser motion via the full nonlinear model

error = S_c-S_rc;                          % Error via the Encke method
error_n = S_c-S_rcn;                       % Error via the full nonlinear model

e(:,1) = sqrt(dot(error, error, 2));       % State error (L2 norm) via Encke's method
e(:,2) = sqrt(dot(error_n, error_n, 2));   % State error (L2 norm) via Newtonian formulation

save NumericalError s0 S Sn;

%% Plotting and results %% 
% Plot results 
figure
view(3) 
hold on
plot3(S_c(:,1), S_c(:,2), S_c(:,3), 'y', 'Linewidth', 0.9); 
plot3(S_rcn(:,1), S_rcn(:,2), S_rcn(:,3), 'r', 'Linewidth', 0.9); 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'b', 'Linewidth', 0.9); 
hold off
legend('True trajectory', 'Newton', 'Encke', 'Location', 'northeast'); 
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;

figure 
hold on
plot(t, log(e(:,1)), 'b'); 
plot(t, log(e(:,2)), 'r');
hold off
grid on
xlabel('$t$'); 
ylabel('Absolute error $\log{e}$');
legend('Encke', 'Newton');
