%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 23/03/21 % 

%% Error study %% 
% This script provides an interface to the evolution of the relative dynamics for differente initial 
% conditions. 

% The relative motion of two spacecraft in two halo orbits around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion % 

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

% Relative motion happens between a halo orbit for the target and random
% initial conditions (pseudo-halo) for the chaser. 

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
%Time span 
dt = 1e-3;                          %Time step
tmax = pi;                          %Maximum time of integration (corresponding to a synodic period)
tspan = 0:dt:tmax;                  %Integration time span

%CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
L = libration_points(mu);           %System libration points
Lem = 384400e3;                     %Mean distance from the Earth to the Moon
T = 2.361e6;                        %Mean period of the Moon around the Earth

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
%Initial relative state 
l = dimensionalizer(Lem, T, Lem/T, 1, 'Position', 0);           %Normalize 1 m in the synodic frame
v = dimensionalizer(Lem, T, Lem/T, 1, 'Velocity', 0);           %Normalize 1 m/s in the synodic frame
e = [l l l v v v];                                              %Initial error
r_t0 = target_orbit.Trajectory(1,1:6);                          %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6)+e;                        %Initial chaser conditions 
rho0 = r_c0-r_t0;                                               %Initial relative conditions
s0 = [r_t0 rho0];                                               %Initial conditions of the target and the relative state

%Integration of the model
[~, S_c] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, r_c0, options);
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%Error in the computation
error = S_c-S_rc;                                           %Position error (L2 norm) via Encke's method          
ep = zeros(size(error,1), 1);                               
for i = 1:size(error,1)
    ep(i) = norm(error(i,:));
end

%% Results %%
%Plot results
figure(1) 
view(3) 
hold on
plot3(S(:,7), S(:,8), S(:,9));  
hold off 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative orbit');

figure(2) 
hold on
plot(t, S(:,7));  
plot(t, S(:,8));  
plot(t, S(:,9));  
hold off 
legend('Synodic x coordinate', 'Synodic y coordinate', 'Synodic z coordinate');
grid on;
title('Relative synodic orbit evolution');

figure(3) 
hold on
plot(t, S(:,10));  
plot(t, S(:,11));  
plot(t, S(:,12));  
hold off 
legend('Synodic x velocity', 'Synodic y velocity', 'Synodic z velocity');
grid on;
title('Relative synodic velocity evolution');

figure(4) 
hold on
plot(t, log(ep), 'b'); 
hold off
grid on
xlabel('Nondimensional epoch'); 
ylabel('Relative error in the synodic frame');
title('Error in the phase space vector (L2 norm)')
legend('Encke method error');