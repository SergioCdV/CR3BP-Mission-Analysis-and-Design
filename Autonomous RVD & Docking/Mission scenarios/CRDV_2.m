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

%CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
Lem = 384400e3;                     %Mean distance from the Earth to the Moon

%% Initial conditions and halo orbit computation %%
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
num = 2;                                                    %Number of orbits to continuate
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

[target_orbit, state_energy] = continuation(num, method, algorithm, object, corrector, setup);

%Generate the NRHO
s0 = [target_orbit.Seeds(end,:).'; reshape(eye(6), [36 1])];
tspan = 0:dt:target_orbit.Period(end);
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);
Sn = S; 

%% Compute the tranfer trajectory
%Parking orbit definition 
R = [1-mu; 0; 0];                       %Primary location in the configuration space
branch = 'L';                           %Manifold branch to globalize
map = 'Secondary primary';              %Poincaré map to use
event = @(t,s)sp_crossing(t,s,mu);      %Integration event

hd = dimensionalizer(Lem, 1, 1, 2000e3, 'Position', 0);                  %Parking orbit altitude

%Integrate the stable manifold backwards and check if it intersects the whereabouts of the parking orbit
manifold = 'S';                                                          %Integrate the stable manifold
seed = S;                                                                %Periodic orbit seed
tspan = 0:1e-3:target_orbit.Period;                                      %Original integration time
rho = 50;                                                                %Density of fibres to analyze
S = invariant_manifold(mu, manifold, branch, seed, rho, tspan, map);     %Initial trajectories

%Relative distance to the primary of interest
distance = zeros(rho,1);    
for i = 1:rho
    %Distance to the orbital altitude
    distance(i) = norm(shiftdim(S.Trajectory(i,S.ArcLength(i),1:3))-R)-hd;  
end

[~, index] = sort(distance);                            %Select the closest manifold to the parking orbit
s0 = shiftdim(S.Trajectory(index(1),1,:));              %Initial conditions to correct
Phi = eye(n);                                           %Initial STM 
Phi = reshape(Phi, [n^2 1]);                            %Initial STM 
sHalo = seed(S.Index(index(1)),1:n).';                  %Halo insertion point
s0 = [s0(1:3); s0(4:6); Phi];                           %Initial conditions
TOF = S.TOF(index(1));                                  %Time of flight
dt = 1e-3;                                              %Time step
tspan = TOF:-dt:0;                                                                  %Integration time
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', event);             %Integration tolerances  
[~, St] = ode113(@(t,s)cr3bp_equations(mu, 1, true, t, s), tspan, s0, options);     %Natural trajectory

%% Modelling in the synodic frame %%
r_t0 = Sn(mod(size(St,1), size(Sn,1)), 1:6);                %Initial target conditions
r_c0 = St(1,1:6);                                           %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Insertion plot 
figure(1)
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3));
plot3(St(:,1), St(:,2), St(:,3));
scatter3(r_c0(1), r_c0(2), r_c0(3))
scatter3(r_t0(1), r_t0(2), r_t0(3))
hold off
legend('Target orbit', 'Target at HOI', 'Chaser at HOI'); 
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
title('Reconstruction of the natural chaser motion');

%% First phase: long-range rendezvous using the TITA approach
%Time of flight 
tf(1) = 0.2;                                            %Nondimensional maneuver end time 
sd = [dimensionalizer(Lem, 1, 1, 1e3, 'Position', 0) 0 0];     %Desired state of the system

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
[St1, dV, state] = TITA_control(mu, tf(1), s0, tol, cost_function, sd, two_impulsive, penalties, target_points, thruster_model);
  
%Long range plot 
figure(1)
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3));
plot3(St1(:,1)+St1(:,7), St1(:,2)+St1(:,8), St1(:,3)+St1(:,9));
hold off
legend('Target orbit', 'TITM arc'); 
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
title('Reconstruction of the natural chaser motion');

%% Second phase: close-range rendezvous
%Phase definition 
tf(2) = ; 

%GNC algorithms definition 
GNC.Algorithms.Guidance = '';               %Guidance algorithm
GNC.Algorithms.Navigation = '';             %Navigation algorithm
GNC.Algorithms.Control = 'SMC';             %Control algorithm
GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                  %Dimension of the control law
GNC.System.mu = mu;                         %System reduced gravitational parameter

%Controller parameters
GNC.Control.SMC.Parameters = SMC_optimization(); 

%Integration time 
tspan = 0:dt:tf(2)-tf(1); 

%Initial conditions 
s0 = St1(end,:); 

%Re-integrate trajectory
[~, St2] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan, s0, options);

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
