%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 07/05/21
% File: HOI_transfer.m 
% Issue: 0 
% Validated: 

%% Halo Orbit Insertion transfer %% 
% This script provides an interface to test the differential corrector aimed to generate HOI trajectories.

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
%Phase space dimension 
n = 6; 

%Time span 
dt = 1e-3;                          %Time step
tf = 2*pi;                          %Rendezvous time
tspan = 0:dt:tf;                    %Integration time spans

%CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
L = libration_points(mu);           %System libration points
Lem = 384400e3;                     %Mean distance from the Earth to the Moon

%Differential corrector set up
maxIter = 20;                       %Maximum number of iterations
tol = 1e-10;                        %Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
%Halo characteristics 
Az = 200e6;                                                         %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Generation of the transfer trajectory
%Definition of the parking orbit 
parking_orbit.Primary = 'First';                                                %Parking orbit primary
parking_orbit.Altitude = dimensionalizer(Lem, 1, 1, 36000e3, 'Position', 0);    %Parking orbit altitude
parking_orbit.Theta = -pi/2;                                                    %Parking orbit true anomaly at TTI

%Redefinition with addtional parameters of the target orbit 
target_orbit.tspan = tspan;                                                     %Original integration time

%Transfer orbit
[transfer_orbit, ~] = transfer_correction('HOI transfer', mu, parking_orbit, target_orbit, maxIter, tol);

%% Plot results 
figure(1) 
view(3)
hold on 
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3));
plot3(transfer_orbit.Trajectory(:,1), transfer_orbit.Trajectory(:,2), transfer_orbit.Trajectory(:,3));
hold off
grid on;
title('Halo Orbit Insertion transfer orbit')
legend('Halo target orbit', 'Parking orbit')
