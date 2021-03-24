%% Phase space study %% 
% This scripts provides an interface to test how the phase space volumens evolve. 

% The relative motion of two spacecraft in two halo orbits around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion 

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
tmax = 4*pi;                        %Maximum time of integration (corresponding to a synodic period)
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

%Continuate the firts halo to locate the chaser spacecraft
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                    %Number of orbits to continuate
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};         %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
direction = 1;                                              %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                         %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(2,1:6), maxIter, tol);

%% Modelling in the synodic frame %% 
r_t0 = target_orbit.Trajectory(1,1:6);                      %Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0];                                           %Initial conditions of the target and the relative state
ds0 = [1 0 0 0 0 0];                                        %Initial phase space error vector
s0 = [s0 ds0];                                              %Initial conditions of the target and the relative state and the error

%Integration of the model
[~, S_c] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, r_c0, options);
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke V', t, s), tspan, s0, options);

%% Results %% 
% Plot results 
figure(1) 
view(3) 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
hold on
for i = 1:size(S_rc,1) 
end
hold off
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Phase space evolution');
