%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 25/04/21 % 

%% GNC 11: Complete rendezvous mission example 1 %% 
% This script provides an interface to test the general control scheme for a rendezvous, docking and undocking mission. 

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
tf = [0.6 pi 2*pi];                 %Switching times
tspan = 0:dt:tf(3);                 %Integration time span

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
Az = 200e6;                                                 %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Modelling in the synodic frame %%
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% First phase: long-range rendezvous using the TITA approach
%Differential corrector set up
tol = 1e-5;                                 %Differential corrector tolerance

%Cost function matrices
penalties.R = eye(3);                       %Penalty on the impulse
penalties.Q = eye(6);                       %Penalty on the state error
penalties.M  = 0.1*eye(6);                  %Penalty on the state noise

%Select measuring times 
target_points.Noise = true;                 %Boolean to account for state noise
target_points.Times = tf(1)*rand(1,3);      %Times to measure the state noise

thruster_model.Sigma = 0.01;                %Velocity noise dependance on the velocity impulse
thruster_model.Rotation = eye(3);           %Rotational misalignment of the thrusters

%Cost function 
cost_function = 'Position';                 %Cost function to target
two_impulsive = true;                       %Two-impulsive rendezvous boolean

%TITA controller
[St1, dV, state] = TITA_control(mu, tf(1), s0, tol, cost_function, zeros(1,3), two_impulsive, ...
                                penalties, target_points, thruster_model);

%% Second phase: docking and coordinated flight 
%GNC algorithms definition 
GNC.Algorithms.Guidance = '';               %Guidance algorithm
GNC.Algorithms.Navigation = '';             %Navigation algorithm
GNC.Algorithms.Control = 'SMC';             %Control algorithm
GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                  %Dimension of the control law
GNC.System.mu = mu;                         %System reduced gravitational parameter
GNC.Control.SMC.Parameters = [1 1 0.9 0.1]; %Controller parameters

%Integration time 
tspan = 0:dt:tf(2)-tf(1); 

%Initial conditions 
s0 = St1(end,:); 

%Re-integrate trajectory
[~, St2] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan, s0, options);

%% Third phase: undocking 
constraint.Constrained = false;         %No constraints on the maneuver
constraint.SafeDistance = 1e-4;         %Safety distance at the collision time

[St3, dV3, tm] = FMSC_control(mu, tf(3)-(tf(2)+dt), Sn(end,1:6), St2(end,1:12), eye(3), tol, constraint, 'Best');

%% Final results 
St = [St1(:,1:12); St2(:,1:12); St3];   %Complete trajectory 
tspan = 0:dt:tf(3)+dt;                  %Total integration time

%% Plotting
figure(2) 
view(3) 
hold on 
plot3(St(:,7), St(:,8), St(:,9), 'r'); 
hold off
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative motion trajectory in the configuration space');

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan(1:size(St,1)), St(:,7)); 
plot(tspan(1:size(St,1)), St(:,8)); 
plot(tspan(1:size(St,1)), St(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('x coordinate', 'y coordinate', 'z coordinate');
title('Relative position evolution');
subplot(1,2,2)
hold on
plot(tspan(1:size(St,1)), St(:,10)); 
plot(tspan(1:size(St,1)), St(:,11)); 
plot(tspan(1:size(St,1)), St(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinate');
grid on;
legend('x velocity', 'y velocity', 'z velocity');
title('Relative velocity evolution');
