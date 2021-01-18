%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 09/01/20
% File: test_2.m 
% Issue: 0 
% Validated: 

%% Test 2 %%
% This scripts provides a test interface for the rest of the library
% functions. Test 2 is concerned with differential correction algorithms validation.

%% Test values and constants
%Initial conditions
mu = 0.0121505856;                                          %Reduced gravitational parameter of the system (Earth-Moon)
L = libration_points(mu);                                   %System libration points
Az = 200e6;                                                 %Orbit amplitude out of the synodic plane
Ax = 200e6;                                                 %Orbit amplitude in the synodic plane 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ax = dimensionalizer(384400e3, 1, 1, Ax, 'Position', 0);    %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around Li
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param1 = [-1 Az Ln gamma m];                                %Halo orbit parameters
param2 = [Ax Az 0 pi/2 Ln gamma m];                         %Lyapunov orbit parameters

%Correction parameters 
maxIter = 50;      %Maximum allowed iterations
tol = 1e-5;       %Tolerance 

%% Functions
%Set graphical environment 
set_graphics(); 

%Compute seeds
[halo_seed, haloT] = object_seed(mu, param1, 'Halo');
lyapunov_seed = object_seed(mu, param2, 'Lyapunov');
axial_seed = [1.1389 0 0 0 -0.2647 0.3659];
vertical_seed = [1.0613 0 0 0 -1.9891 0.4740];
butterfly_seed = [1.0406 0 0.1735 0 -0.0770 0];

%Lyapunov orbit
[lyapunov_orbit, state(1)] = differential_correction('Planar', mu, lyapunov_seed, maxIter, tol);

%Halo orbit
[halo_orbit, state(2)] = differential_correction('Periodic MS', mu, halo_seed, maxIter, tol, 15, haloT);

%Distant Retrograde Orbit 
[dro_orbit, state(3)] = differential_correction('Double Symmetric', mu, halo_seed, maxIter, tol);

%Axial Orbit
[axial_orbit, state(4)] = differential_correction('Axis Symmetric', mu, axial_seed, maxIter, tol);

%Vertical Orbit 
[vertical_orbit, state(5)] = differential_correction('Axis Symmetric', mu, vertical_seed, maxIter, tol);

%Butterfly Orbit
[butterfly_orbit, state(6)] = differential_correction('Plane Symmetric', mu, butterfly_seed, maxIter, tol);

%% Plotting
figure(1) 
plot3(lyapunov_seed(:,1), lyapunov_seed(:,2), lyapunov_seed(:,3));
xlabel('Synodic normalized x coordinate');
ylabel('Synodic normalized y coordinate');
zlabel('Synodic normalized z coordinate');
title('Lyapunov orbit seed');
grid on;

figure(2) 
plot3(lyapunov_orbit.Trajectory(:,1), lyapunov_orbit.Trajectory(:,2), lyapunov_orbit.Trajectory(:,3));
xlabel('Synodic normalized x coordinate');
ylabel('Synodic normalized y coordinate');
zlabel('Synodic normalized z coordinate');
title('Converged Lyapunov orbit');
grid on;

figure(3) 
plot3(halo_orbit.Trajectory(:,1), halo_orbit.Trajectory(:,2), halo_orbit.Trajectory(:,3));
xlabel('Synodic normalized x coordinate');
ylabel('Synodic normalized y coordinate');
zlabel('Synodic normalized z coordinate');
title('Halo orbit');
grid on;

figure(4)
plot3(dro_orbit.Trajectory(:,1), dro_orbit.Trajectory(:,2), dro_orbit.Trajectory(:,3));
xlabel('Synodic normalized x coordinate');
ylabel('Synodic normalized y coordinate');
zlabel('Synodic normalized z coordinate');
title('Converged DRO');
grid on;

figure(5) 
plot3(axial_orbit.Trajectory(:,1), axial_orbit.Trajectory(:,2), axial_orbit.Trajectory(:,3));
xlabel('Synodic normalized x coordinate');
ylabel('Synodic normalized y coordinate');
zlabel('Synodic normalized z coordinate');
title('Converged axial orbit');
grid on;

figure(6) 
plot3(vertical_orbit.Trajectory(:,1), vertical_orbit.Trajectory(:,2), vertical_orbit.Trajectory(:,3));
xlabel('Synodic normalized x coordinate');
ylabel('Synodic normalized y coordinate');
zlabel('Synodic normalized z coordinate');
title('Converged vertical orbit');
grid on;

figure(7) 
plot3(butterfly_orbit.Trajectory(:,1), butterfly_orbit.Trajectory(:,2), butterfly_orbit.Trajectory(:,3));
xlabel('Synodic normalized x coordinate');
ylabel('Synodic normalized y coordinate');
zlabel('Synodic normalized z coordinate');
title('Converged butterfly orbit');
grid on;

%% Auxiliary functions 
function set_graphics()
    %Set graphical properties
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
    set(groot, 'defaultAxesFontSize', 11); 
    set(groot, 'defaultAxesGridAlpha', 0.3); 
    set(groot, 'defaultAxesLineWidth', 0.75);
    set(groot, 'defaultAxesXMinorTick', 'on');
    set(groot, 'defaultAxesYMinorTick', 'on');
    set(groot, 'defaultFigureRenderer', 'painters');
    set(groot, 'defaultLegendBox', 'off');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultLegendLocation', 'best');
    set(groot, 'defaultLineLineWidth', 1); 
    set(groot, 'defaultLineMarkerSize', 3);
    set(groot, 'defaultTextInterpreter','latex');
end