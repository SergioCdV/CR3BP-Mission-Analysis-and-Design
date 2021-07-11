%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 14/02/21
% File: test_5.m 
% Issue: 0 
% Validated: 

%% Test 5 %%
% This scripts provides a test interface for the rest of the library
% functions. 

% Test 5 is concerned to the globalization of NHROs and their unstable and stable manifolds.

%% Test values and constants
%Set graphical environment 
set_graphics(); 

%Numerical setup 
format long; 
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Integration tolerances

%Initial conditions
mu = 0.0121505856;                                          %Reduced gravitational parameter of the system (Earth-Moon)
L = libration_points(mu);                                   %System libration points
Az = 195e6;                                                 %Orbit amplitude out of the synodic plane. Play with it! 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param = [-1 Az Ln gamma m];                                 %Halo orbit parameters (-1 being for southern halo)

%Correction parameters 
dt = 1e-3;                                                  %Time step to integrate converged trajectories
maxIter = 20;                                               %Maximum allowed iterations in the differential correction schemes
tol = 1e-10;                                                %Differential correction tolerance 
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                   %Number of orbits to continuate
direction = -1;                                             %Direction to continuate (to the Earth)
   
%% Functions
%Compute the NRHO
[halo_seed, haloT] = object_seed(mu, param, 'Halo');        %Generate a halo orbit seed

%Continuation procedure 
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, haloT};                       %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
setup = [mu maxIter tol direction];                         %General setup

[Results_energy, state_energy] = continuation(num, method, algorithm, object, corrector, setup);

%% Generate the NRHO
s0 = [Results_energy.Seeds(end,:).'; reshape(eye(6), [36 1])];
tspan = 0:dt:Results_energy.Period(end);
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);

%% Globalization of the manifolds
%Manifold definition
Branch = ['L' 'L'];            %Directions to propagate the manifolds
TOF = 2*pi;                    %Time of flight
rho = 50;                      %Manifold fibers to compute 

%Computation flags
long_rendezvous = true;         %Flag to allow for long rendezvous
position_fixed = false;         %Flag to determine a final target state
graphics = true;                %Flag to plot the manifolds

halo_orbit.Trajectory = S;      %Target orbit

%Trajectory design core
[Sg, dV] = HCNC_guidance(mu, Branch, rho, Ln, halo_orbit, TOF, long_rendezvous, position_fixed, graphics);

%% Plotting and results 
%Plot the transfer
figure(1) 
view(3)
hold on 
plot3(halo_orbit.Trajectory(:,1), halo_orbit.Trajectory(:,2), halo_orbit.Trajectory(:,3));
plot3(Sg.Trajectory(:,1), Sg.Trajectory(:,2), Sg.Trajectory(:,3), 'k');
hold off
grid on;
xlabel('Synodic normalized $x$ coordinate');
ylabel('Synodic normalized $y$ coordinate');
zlabel('Synodic normalized $z$ coordinate');
title('Homoclinic rendezvous trajectory') 
legend('Target NHRO', 'Homoclinic trajectory')