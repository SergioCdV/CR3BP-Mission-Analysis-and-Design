%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% Homoclinic connection %% 
% This script provides an interface to test homoclinic rendezvous in the CR3BP.

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

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
tspan = 0:dt:tf;                    %Integration time span
tspann = 0:dt:2*pi;                 %Integration time span

%CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
L = libration_points(mu);           %System libration points
Lem = 384400e3;                     %Mean distance from the Earth to the Moon

%Differential corrector set up
nodes = 10;                         %Number of nodes for the multiple shooting corrector
maxIter = 20;                       %Maximum number of iterations
tol = 1e-10;                        %Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
%Halo characteristics 
Az = 120e6;                                                         %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Homoclinic connection
%Manifold definition
Branch = ['L' 'L'];             %Directions to propagate the manifolds
TOF = tf;                       %Time of flight
rho = 50;                       %Manifold fibers to compute 
 
%Computation flags
long_rendezvous = true;         %Flag to allow for long rendezvous
position_fixed = false;         %Flag to determine a final target state
graphics = true;                %Flag to plot the manifolds

%Final target state
target_orbit.TargetState = shiftdim(target_orbit.Trajectory(1,1:n)); 

%Trajectory design core
[Sg, dV] = HCNC_guidance(mu, Branch, rho, Ln, target_orbit, TOF, long_rendezvous, position_fixed, graphics);

%% Plot results
%Plot the transfer
figure 
view(3)
hold on 
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3));
plot3(Sg.Trajectory(:,1), Sg.Trajectory(:,2), Sg.Trajectory(:,3));
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
text(L(1,Ln)+1e-3, L(2,Ln), 0, '$L_1$');
scatter3(1-mu, L(2,Ln), 0, 'k', 'filled');
text(1-mu+1e-3, L(2,Ln), 0, '$M_2$');
hold off
grid on;
title('Homoclinic connection for an $L_1$ halo orbit')
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
legend('Target halo orbit', 'Homoclinic orbit', 'Location', 'northeast')