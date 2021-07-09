%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 20/02/21 % 

%% Linear relative motion models %% 
% This script provides a complete script dedicated to linear models 
% for relative motion in the CR3BP.

% The relative motion of two spacecraft in a halo orbit around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion % 

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
%Time span 
dt = 1e-3;                          %Time step
tmax = 0.8*pi;                      %Maximum time of integration (corresponding to a synodic period)
tspan = 0:dt:tmax;                  %Integration time span

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
Az = 120e6;                                                 %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
cn = legendre_coefficients(mu, Ln, gamma, 2);               %Legendre coefficients of interest                                   
m = 1;                                                      %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[halo_orbit, state] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Modelling in the synodic frame %% 
r_t0 = halo_orbit.Trajectory(1,1:6);                        %Initial target conditions
r_c0 = halo_orbit.Trajectory(2,1:6);                        %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0];                                           %Initial conditions of the target and the relative state

%Integration of the models
[~, S_c] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, r_c0, options);

[~, S_t] = ode113(@(t,s)lr_model(mu, cn, true, false, 'Target', t, s), tspan, s0, options);
[~, S_fl] = ode113(@(t,s)lr_model(mu, cn, true, false, 'Fixed libration', t, s), tspan, s0, options);
[t, S_ml] = ode113(@(t,s)lr_model(mu, cn, true, false, 'Moving libration', t, s), tspan, s0, options);

%Reconstructed chaser motion 
S_ct = S_t(:,1:6)+S_t(:,7:12);                              %Chaser motion from the target linearization
S_cfl = S_fl(:,1:6)+S_fl(:,7:12);                           %Chaser motion from the fixed libration point linearization
S_cml = S_ml(:,1:6)+S_ml(:,7:12);                           %Chaser motion from the moving libration point linearization

%% Error in the approximations %%
%Preallocation 
error_t = zeros(1,size(S_c,1));                             %Error using the target linearization model
error_fl = zeros(1,size(S_c,1));                            %Error using the fixed libration point linearization model
error_ml = zeros(1,size(S_c,1));                            %Error using the moving libration point linearization model

%Main computation
for i = 1:size(S_c,1)
    error_t(i) =  norm(S_c(i,:)-S_ct(i,:));                 %Error using the target linearization model
    error_fl(i) = norm(S_c(i,:)-S_cfl(i,:));                %Error using the fixed libration point linearization model
    error_ml(i) = norm(S_c(i,:)-S_cml(i,:));                %Error using the moving libration point linearization model
end

%% Results %% 
%Plot results 
figure(1) 
view(3) 
hold on
plot3(S_c(:,1), S_c(:,2), S_c(:,3), 'm'); 
plot3(S_ct(:,1), S_ct(:,2), S_ct(:,3), 'k'); 
plot3(S_cfl(:,1), S_cfl(:,2), S_cfl(:,3), 'b'); 
plot3(S_cml(:,1), S_cml(:,2), S_cml(:,3), 'r'); 
hold off
legend('Expected trajectory', 'RLM', 'SRLLM', 'URLLM', 'Location', 'northeast'); 
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $x$ coordinate');
zlabel('Synodic $x$ coordinate');
grid on;
title('Reconstruction of the chaser motion via the proposed linear models');

figure(2) 
hold on
plot(t, log(error_t), 'k'); 
plot(t, log(error_ml), 'r'); 
plot(t, log(error_fl), 'b'); 
hold off
grid on
xlabel('Nondimensional epoch'); 
ylabel('Absolute error $\log{e}$');
legend('RLM', 'SRLLM', 'URLLM');
title('Absolute error between the true and linear dynamics')


