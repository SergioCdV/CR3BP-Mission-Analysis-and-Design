%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/06/21
% File: test_6.m 
% Issue: 0 
% Validated: 

%% Test 6 %%
% This scripts provides a test interface for the rest of the library
% functions. 

% Test 6 is concerned with the computation of libration curves and
% subsitutes.

%% Test values and constants
%Set graphical environment 
set_graphics(); 

%Numerical setup 
format long; 
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Integration tolerances

%Initial conditions
mu = 0.0121505856;                                          %Reduced gravitational parameter of the system (Earth-Moon)
L = libration_points(mu);                                   %System libration points
Az = 10e6;                                                  %Orbit amplitude out of the synodic plane. Play with it! 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param = [-1 Az Ln gamma m];                                 %Halo orbit parameters (-1 being for southern halo)
dt = 1e-3;                                                  %Time step to integrate converged trajectories

%% Functions
%Compute the NRHO
[halo_seed, haloT] = object_seed(mu, param, 'Halo');        %Generate a halo orbit seed

%% Libration curves
Lc = libration_curves(mu, halo_seed(:,1:3), Ln); 

%% Plotting and results 
%Plot the transfer
figure 
view(3)
hold on 
plot3(Lc(:,1), Lc(:,2), Lc(:,3));
hold off
grid on;
title('Libration curves')