%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 20/02/21 % 

%% First contact %% 
% This scripts provides a first contact with the RVD project. The objective
% is to qualitatively analyze whether the Clohessy-Whiltshire
% relative motion linear model could be use in the CR3BP problem. 

% The relative motion of two spacecraft in a halo orbit around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and with the CW model % 

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frames as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode45)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
% Time span 
dt = 1e-3;                  %Time step
tmax = 2*pi;                %Maximum time of integration (corresponding to a synodic period)
tspan = 0:dt:tmax;          %Integration time span

% CR3BP constants 
mu = 0.0121505;             %Earth-Moon reduced gravitational parameter
L = libration_points(mu);   %System libration points
Lem = 384400e3;             %Mean distance from the Earth to the Moon

% CR3BP integration flags 
flagVar = 1;                %Integrate the dynamics and the first variotional equations 
direction = 1;              %Integrate forward in time 

% Differential corrector set up
nodes = 10;                 %Number of nodes for the multiple shooting corrector
maxIter = 20;               %Maximum number of iterations
tol = 1e-7;                 %Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
%Halo characteristics 
Az = 75e6;                                                  %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[halo_orbit, state] = differential_correction('Jacobi Constant MS', mu, halo_seed, maxIter, tol, nodes, period, 0);

%% Modelling %% 

%% Results %% 