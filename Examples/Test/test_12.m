%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 04/07/22
% File: test_12.m 
% Issue: 0 
% Validated: 

%% Test 12 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 12 is concerned with differential correction algorithms validation and the exploration of 
% different families of relative orbits

% Credit to Grebow, 2006, for his initial seeds!

set_graphics();          % Set graphical environment 

%% Test values and constants
% Initial conditions
mu = 0.0121505856;                                          % Reduced gravitational parameter of the system (Earth-Moon)
Lem = 384400e3;                                             % Characteristic distance of the system (Earth-Moon)
L = libration_points(mu);                                   % System libration points
Az = 20e5;                                                  % Orbit amplitude out of the synodic plane. Play with it!
Ax = 20e6;                                                  % Orbit amplitude in the synodic plane. Play with it! 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         % Normalize distances for the E-M system
Ax = dimensionalizer(Lem, 1, 1, Ax, 'Position', 0);         % Normalize distances for the E-M system
Ln = 2;                                                     % Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          % Li distance to the second primary
m = 1;                                                      % Number of periods to compute
param_halo = [1 Az Ln gamma m];                             % Halo orbit parameters (-1 for southern halo)
param_lyap = [Ax Az 0 0 Ln gamma m];                        %Lyapunov orbit parameters

% Correction parameters 
maxIter = 50;                                               % Maximum allowed iterations in the differential correction schemes
tol = 1e-10;                                                % Tolerance 

%% Functions
% Compute absolute seeds
[halo_seed, haloT] = object_seed(mu, param_halo, 'Halo');       % Generate a halo orbit seed
lyapunov_seed = object_seed(mu, param_lyap, 'Lyapunov');        % Generate a Lyapunov orbit seed
axial_seed = [0.8431 0 0 0 0.1874 0.4000];                      % State vector of an axial orbit
vertical_seed = [0.9261 0 0.3616 0 -0.0544  0];                 % State vector of a vertical orbit
butterfly_seed = [1.0406 0 0.1735 0 -0.0770 0];                 % State vector of a butterfly orbit

% Lyapunov orbit
[lyapunov_orbit, state(1)] = differential_correction('Planar', mu, lyapunov_seed, maxIter, tol);

param_lyap = [Ax Az 0 0 Ln gamma 20];                           % Lyapunov orbit parameters
lyapunov_seed = object_seed(mu, param_lyap, 'Lyapunov');        % Generate a Lyapunov orbit seed

% Halo orbit 
[halo_orbit, state(2)] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

% Distant Retrograde Orbit (only for L2)
[dro_orbit, state(3)] = differential_correction('Planar', mu, halo_seed, maxIter, tol);

% Axial Orbit
[axial_orbit, state(4)] = differential_correction('Axis Symmetric', mu, axial_seed, maxIter, tol);

% Vertical Orbit 
[vertical_orbit, state(5)] = differential_correction('Double Plane Symmetric', mu, vertical_seed, maxIter, tol);

% Butterfly Orbit
[butterfly_orbit, state(6)] = differential_correction('Plane Symmetric', mu, butterfly_seed, maxIter, tol);

%% Computation of the relative orbit
% Continuation and chaser trajectory 
Bif_tol = 1e-2;                                             % Bifucartion tolerance on the stability index
num = 5;                                                    % Number of orbits to continuate
direction = 1;                                              % Direction to continuate (to the Earth)
method = 'SPC';                                             % Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                % Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, halo_orbit.Period};           % Object and characteristics to continuate
corrector = 'Plane Symmetric';                              % Differential corrector method
setup = [mu maxIter tol direction];                         % General setup

[Results_energy, ~] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, Results_energy.Seeds(end,:), maxIter, tol);

Sc = chaser_orbit.Trajectory;                               % Chaser trajectory initial conditions
St = halo_orbit.Trajectory;                                 % Target trajectory initial conditions

% Relative initial conditions 
maxIter = 50;                                               % Maximum number of iterations
s0 = [St(1,1:6) Sc(1,1:6)-St(1,1:6)];                       % Relative initial conditions

% Correction
[rel_orbit, state(7)] = differential_correction('Plane Symmetric', mu, s0(7:12), maxIter, tol);

% New chaser orbit 
Sc = St(:,1:6) + rel_orbit.Trajectory(:,7:12);

%% Plotting
figure(1) 
view(3)
hold on 
plot3(rel_orbit.Trajectory(:,7), rel_orbit.Trajectory(:,8), rel_orbit.Trajectory(:,9), 'b');
hold off
xlabel('Relative synodic $x$ coordinate');
ylabel('Relative synodic $y$ coordinate');
zlabel('Relative synodic $z$ coordinate');
title('Relative periodic orbit');
grid on;

figure(2) 
view(3)
hold on
J = plot3(St(:,1), St(:,2), St(:,3), 'b');
H = plot3(Sc(:,1), Sc(:,2), Sc(:,3), 'r');
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
title('Earth-Moon $L_1$ planar Lyapunov orbit');
legend({'Target orbit', 'Chaser orbit'}, 'Location', 'northeast');
grid on;