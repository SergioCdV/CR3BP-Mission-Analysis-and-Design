%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 20/02/21 % 

%% First contact %% 
% This scripts provides a first contact with the RVD project. 

% The relative motion of two spacecraft in a halo orbit around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion % 

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode45)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
% Time span 
dt = 1e-5;                          %Time step
tmax = 5*pi;                       %Maximum time of integration (corresponding to a synodic period)
tspan = 0:dt:tmax;                  %Integration time span

% CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
L = libration_points(mu);           %System libration points
Lem = 384400e3;                     %Mean distance from the Earth to the Moon

% CR3BP integration flags 
flagVar = 1;                        %Integrate the dynamics and the first variotional equations 
direction = 1;                      %Integrate forward in time 

% Differential corrector set up
nodes = 10;                         %Number of nodes for the multiple shooting corrector
maxIter = 20;                       %Maximum number of iterations
tol = 1e-7;                         %Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
%Halo characteristics 
Az = 100e6;                                                 %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
K = 1e-10;                                                   %Nondimensional distance from target to chaser

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[halo_orbit, state] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Modelling in the synodic frame %% 
r_t0 = halo_orbit.Trajectory(1,1:6);                        %Initial target conditions
r_c0(1:3) = r_t0(1:3)+K*rand(1,3);                          %Initial chaser conditions 
rho0 = [K*rand(1,3) r_t0(4:6)];                             %Initial relative conditions
s0 = [r_t0 rho0];                                           %Initial conditions of the target and the relative state

%Integration of the model
[~, S_c] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, r_c0, options);
[t, S] = ode113(@(t,s)fulrel_motion(mu, true, false, true, t, s), tspan, s0, options);
%[t, Sn] = ode113(@(t,s)fulrel_motion(mu, true, false, false, t, s), tspan, s0, options);

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method
error = S_c-S_rc;                                           %Compute the error via the full nonlinear model with the Encke method
%S_rcn = Sn(:,1:6)+Sn(:,7:12);                               %Reconstructed chaser motion via the full nonlinear model
%error_n = S_c-S_rce;                                        %Compute the error via the full nonlinear model

%% Results in the inertial frame %% 
%Preallocation 
inertial_error = zeros(size(S_c,1),3);      %Position error in the inertial frame

%Main computation
for i = 1:size(S_c,1)
     elaps = dt*(i-1);                                                       %Elaps time
     inertial_error(i,:) = (inertial2synodic(elaps, error(i,1:3).', 1)).';   %Position error in the synodic frame
end

%% Results in the libration synodic frame %% 
%Preallocation 
libration_error = zeros(size(S_c,1),3);     %Position error in the synodic frame

%Main computation
for i = 1:size(S_c, 1)
     libration_error(i,:) = (synodic2lagrange(mu, gamma, Ln, error(i,1:3).', 1)).'; %Position error in the synodic frame
end

%% Results in the Frenet-Serret frame of the target orbit %% 
%Preallocation
T = zeros(size(S_c,1),3,3);             %Synodic to Frenet frame rotation matrix
frenet_error = zeros(size(S_c,1), 3);   %Position error in the Frenet frame

%Main computation
for i = 1:size(S_c,1)
    T(i,1:3,1:3) = frenet_triad(mu, S(i,1:6));                   %Frenet-Serret frame
    frenet_error(i,:) = (shiftdim(T(i,:,:)).'*error(i,1:3).').'; %Position error in the Frenet-Serret frame of the target orbit
end

%% Results %% 
% Plot results 
figure(1) 
view(3) 
hold on
plot3(S(:,1), S(:,2), S(:,3)); 
plot3(S_c(:,1), S_c(:,2), S_c(:,3)); 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
hold off
legend('Target motion', 'Chaser motion', 'New chaser motion'); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;

figure(2) 
hold on
plot(t, error); 
%plot(t, error_n(:,1), 'r');
hold off
grid on
xlabel('Nondimensional epoch'); 
ylabel('Relative error in the synodic frame');
title('Error in the chaser motion reconstruction')

figure(3) 
plot3(inertial_error(:,1), inertial_error(:,2), inertial_error(:,3)); 
xlabel('Nondimensional x coordinate'); 
ylabel('Nondimensional y coordinate');
zlabel('Nondimensional z coordinate');
grid on
title('Position error in the inertial frame'); 

figure(4) 
plot3(libration_error(:,1), libration_error(:,2), libration_error(:,3)); 
xlabel('Nondimensional x coordinate'); 
ylabel('Nondimensional y coordinate');
zlabel('Nondimensional z coordinate');
grid on
title('Position error in the libration synodic frame'); 

figure(5) 
plot3(frenet_error(:,1), frenet_error(:,2), frenet_error(:,3)); 
xlabel('Nondimensional x coordinate'); 
ylabel('Nondimensional y coordinate');
zlabel('Nondimensional z coordinate');
grid on
title('Position error in the Frenet-Serret frame');
