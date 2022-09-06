%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/02/21
% File: test_3.m 
% Issue: 0 
% Validated: 

%% Test 3 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 3 is concerned with the computation of invariant manifols 
% associated with periodic orbits

%% General setup 
format long;                                                % Set command window numerical format
set_graphics();                                             % Set up the graphical environment

%% Initial conditions 
mu = 0.0121505856;                                          % Reduced gravitational parameter of the system
L = libration_points(mu);                                   % Libration points computation

%Generate a Lyapunov orbit seed
mu = 0.0121505856;                                          % Reduced gravitational parameter of the system (Earth-Moon)
Lem = 384400e3;                                             % Characteristic distance of the system (Earth-Moon)
L = libration_points(mu);                                   % System libration points
Az = 20e5;                                                  % Orbit amplitude out of the synodic plane. Play with it!
Ax = 20e6;                                                  % Orbit amplitude in the synodic plane. Play with it! 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         % Normalize distances for the E-M system
Ax = dimensionalizer(Lem, 1, 1, Ax, 'Position', 0);         % Normalize distances for the E-M system
Ln = 2;                                                     % Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          % Li distance to the second primary
m = 1;                                                      % Number of periods to compute
param_halo = [1 Az Ln gamma m];                             % Halo orbit parameters (-1 for southern halo)
param = [Ax Az 0 0 Ln gamma m];                             % Lyapunov orbit parameters
lyapunov_seed = object_seed(mu, param, 'Lyapunov');         % Lyapunov seed     
[halo_seed, haloT] = object_seed(mu, param_halo, 'Halo');   % Halo orbit seed

% Time integration
Tf = pi;                                                    % Final time
dt = 1e-3;                                                  % Time step
tspan = 0:dt:Tf;                                            % Integration time span

% Differential correction scheme set up
maxIter = 20;                                               % Maximum allowed iterations in the differential correction schemes
tol = 1e-5;                                                 % Tolerance 

%% Main test
% Orbit computation
[lyapunov_orbit, state] = differential_correction('Planar', mu, lyapunov_seed, maxIter, tol);
[halo_orbit, state(2)] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

% Manifolds computation
rho = 10;                    % Number of manifold fibers to compute

manifold_ID = 'S';           % Stable manifold (U or S)
manifold_branch = 'L';       % Left branch of the manifold (L or R)

StableManifold = invariant_manifold(mu, Ln, manifold_ID, manifold_branch, halo_orbit.Trajectory, rho, tspan);

manifold_ID = 'U';           % Unstable manifold (U or S)
manifold_branch = 'L';       % Left branch of the manifold (L or R)

UnstableManifold = invariant_manifold(mu, Ln, manifold_ID, manifold_branch, halo_orbit.Trajectory, rho, tspan);

%% Plotting and results 
figure(1)
hold on 
plot3(lyapunov_orbit.Trajectory(:,1), lyapunov_orbit.Trajectory(:,2), lyapunov_orbit.Trajectory(:,3), 'k');
plot3(halo_orbit.Trajectory(:,1), halo_orbit.Trajectory(:,2), halo_orbit.Trajectory(:,3), 'k');
for i = 1:size(StableManifold.Trajectory,1)
    ManifoldAux = shiftdim(StableManifold.Trajectory(i,:,:));
    plot3(ManifoldAux(1:StableManifold.ArcLength(i),1), ManifoldAux(1:StableManifold.ArcLength(i),2), ManifoldAux(1:StableManifold.ArcLength(i),3), 'g');
end

for i = 1:size(UnstableManifold.Trajectory,1)
    ManifoldAux = shiftdim(UnstableManifold.Trajectory(i,:,:));
    plot3(ManifoldAux(1:UnstableManifold.ArcLength(i),1), ManifoldAux(1:UnstableManifold.ArcLength(i),2), ManifoldAux(1:UnstableManifold.ArcLength(i),3), 'r');
end
%legend('Lyapunov orbit', 'Halo orbit', 'Halo orbit $W^s$', 'Halo orbit $W^u$', 'AutoUpdate', 'off', 'Location', 'northeast');
scatter(1-mu, 0, 'k', 'filled');
scatter(L(1,Ln), 0, 'k', 'filled');
text([1-mu+1e-2 L(1,Ln)+0.02], [0 0], {'$M_2$', '$L_2$'});
hold off
grid on;
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
title('Unstable and stable manifolds of an $L_2$ Lyapunov orbit');
