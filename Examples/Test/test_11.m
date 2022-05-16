%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 09/04/22
% File: test_11.m 
% Issue: 0 
% Validated: 

%% Test 11 %%
% This scripts provides a test interface for the rest of the library
% functions. 

% Test 11 is concerned with the Encke's method integration for collisions in the CR3BP.

%% Test values and constants
%Set graphical environment 
set_graphics(); 

%Initial conditions
mu = 0.0121505856;                                          %Reduced gravitational parameter of the system (Earth-Moon)
L = libration_points(mu);                                   %System libration points
Az = 60e5;                                                  %Orbit amplitude out of the synodic plane. Play with it!
Ax = 60e6;                                                  %Orbit amplitude in the synodic plane. Play with it! 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ax = dimensionalizer(384400e3, 1, 1, Ax, 'Position', 0);    %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param_halo = [1 Az Ln gamma m];                             %Halo orbit parameters (-1 for southern halo)

%Correction parameters 
maxIter = 50;      %Maximum allowed iterations in the differential correction schemes
tol = 1e-10;       %Tolerance 

%Numerical setup 
format long; 
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Integration tolerances

%% Functions
%Initial conditions and time span 
K = 1;                                   %Number of periods to integrate
dt = 1e-3;                               %Time step
tspan = 0:dt:2*pi;                       %Integration time span

s0 = [1-mu+3e-2; 0; 0; 5e-2; 0; 0];

%Propagation setup 
setup.Method = 'Newton';

%Newton integration
tic
[~, S_N] = ode113(@(t,s)cr3bp_propagator(setup, mu, L, true, false, t, s), tspan, s0, options);
toc

%Encke's integration 
setup.Method = 'Encke';
S0 = s0; 
tic
[t, S_E] = ode113(@(t,s)cr3bp_propagator(setup, mu, L, true, false, t, s), tspan, s0, options);
toc

%% Plotting
%Plot the three orbits
figure(1) 
view(3)
hold on
plot3(S_N(:,1), S_N(:,2), S_N(:,3), 'b')
plot3(S_E(:,1), S_E(:,2), S_E(:,3), 'r')
hold off
xlabel('Synodic normalized $x$ coordinate');
ylabel('Synodic normalized $y$ coordinate');
zlabel('Synodic normalized $z$ coordinate');
title('Numerically integrated halo orbits');
legend('Newton', 'Encke');
grid on;