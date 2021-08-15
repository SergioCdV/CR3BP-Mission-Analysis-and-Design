%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 09/01/21
% File: test_2.m 
% Issue: 0 
% Validated: 

%% Test 9 %%
% This scripts provides a test interface for the rest of the library
% functions. 

% Test 9 is concerned with the symplectic integration of periodic orbits.

% Credit to Grebow, 2006, for his initial seeds!

%% Test values and constants
%Set graphical environment 
set_graphics(); 

%Initial conditions
mu = 0.0121505856;                                          %Reduced gravitational parameter of the system (Earth-Moon)
L = libration_points(mu);                                   %System libration points
Az = 50e6;                                                  %Orbit amplitude out of the synodic plane. Play with it!
Ax = 50e6;                                                  %Orbit amplitude in the synodic plane. Play with it! 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ax = dimensionalizer(384400e3, 1, 1, Ax, 'Position', 0);    %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param_halo = [1 Az Ln gamma m];                             %Halo orbit parameters (-1 for southern halo)

%Correction parameters 
maxIter = 50;      %Maximum allowed iterations in the differential correction schemes
tol = 1e-10;       %Tolerance 

%% Functions
%Compute seeds
[halo_seed, haloT] = object_seed(mu, param_halo, 'Halo');   %Generate a halo orbit seed

%Halo orbit 
[halo_orbit, state(2)] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Symplectic integration 
s0 = halo_orbit.Trajectory(1,1:6);      %Initial conditions 
dt = 1e-3;                              %Time step
tspan = 0:dt:0.1;         %Integration time span 
[S, z] = symplectic_flow(mu, Ln, gamma, s0, tspan);

%% Plotting
figure(1) 
view(3)
hold on
%plot3(halo_orbit.Trajectory(:,1), halo_orbit.Trajectory(:,2), halo_orbit.Trajectory(:,3), 'b');
plot3(S(:,1), S(:,2), S(:,3), 'r')
hold off
xlabel('Synodic normalized $x$ coordinate');
ylabel('Synodic normalized $y$ coordinate');
zlabel('Synodic normalized $z$ coordinate');
title('Halo orbit');
legend('Integrated', 'Symplectic');
grid on;