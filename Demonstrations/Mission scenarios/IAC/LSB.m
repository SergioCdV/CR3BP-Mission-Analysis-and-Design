%% Autonomous RVD and docking in the CR3BP  %%
% Date: 27/08/22

%% Set up
set_graphics(); 
close all

%% Center Manifold Lissajous Shape-based Guidance demonstration for IAC 2022 %% 
% This script provides an interface to demonstrate the LSB guidance core

% Units are non-dimensional and solutions are expressed in the synodic
% reference frame as defined by Howell, 1984.

%% Trajectory generation 
% CR3BP constants 
mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
L = libration_points(mu);           % System libration points
Lem = 384400e3;                     % Mean distance from the Earth to the Moon
T0 = 28*86400/(2*pi);               % Mean period for the Earth-Moon system

% Differential corrector set up
nodes = 10;                         % Number of nodes for the multiple shooting corrector
maxIter = 20;                       % Maximum number of iterations
tol = 1e-10;                        % Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
% Halo characteristics 
Az = 20e6;                                                          % Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 1;                                                             % Orbits around L1
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute a halo orbit seed 
halo_param = [1 Az Ln gamma m];                                     % Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

% Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                                     % Bifucartion tolerance on the stability index
num = 5;                                                            % Number of orbits to continuate
method = 'SPC';                                                     % Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                        % Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};                 % Object and characteristics to continuate
corrector = 'Plane Symmetric';                                      % Differential corrector method
direction = 1;                                                      % Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                                 % General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);

% Halo characteristics 
Az = 10e6;                                                          % Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 1;                                                             % Orbits around L1
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute a halo orbit seed 
halo_param = [1 Az Ln gamma m];                                     % Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a halo orbit seed

[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Setup of the solution method
time_distribution = 'Chebyshev';      % Distribution of time intervals
basis = 'Chebyshev';                  % Polynomial basis to be use
n = [15 15];                          % Polynomial order in the state vector expansion
m = 600;                              % Number of sampling points
cost_function = 'Minimum energy';     % Cost function to be minimized

% Chaser's initial Cartesian state vector
initial_state = chaser_orbit.Trajectory(50,1:6); 

% Target's initial Cartesian state vector
target_state = target_orbit.Trajectory(1000,1:6); 

% Spacecraft propulsion parameters 
T = 5e-4;     % Maximum acceleration 
K = 0;        % Initial input revolutions 

% Setup of the SBOPT routine
options.order = n; 
options.basis = basis;
options.grid = time_distribution; 
options.nodes = m; 
options.cost_function = cost_function;
options.resultsFlag = true; 
options.animations = false;  

%% Results
% Setup of the solution 
GNC.Algorithm = 'SDRE';                 % Solver algorithm
GNC.LQR.StateMatrix = 10*eye(2);        % State error weight matrix
GNC.LQR.ControlMatrix = eye(1);         % Control effort weight matrix
GNC.Tmax = T/sqrt(2)*(T0^2/Lem);        % Constrained acceleration
GNC.TOF = pi;                           % Maneuver time
GNC.SBOPT.setup = options;              % SBOPT setup

% method = 'Prescribed shape-based'; 
% method = 'Dynamics shape-based';
% method = 'Numerical shape-based';
% method = 'Minimum energy';

% Relative solution    
tic
[Sr, u, tf, lissajous_constants] = LSB_guidance(mu, Ln, gamma, initial_state-target_state, method, GNC.TOF, GNC);  
toc 

% Absolute chase trajectory
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
tau = linspace(0,tf,size(Sr,2));
[~, Sc] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tau, target_state, options);
C = Sc.'+Sr;

% Valuation 
[error, merit] = figures_merit(tau, [zeros(size(Sr,2), 6) Sr.']);        % Error performance indices 
effort = control_effort(tau, u, false);                                  % Control effort indices

%% Manifolds computation
rho = 1;                     % Number of manifold fibers to compute
tspan = 0:1e-3:0.1*tf;         % Integration timespan

manifold_ID = 'S';           % Stable manifold (U or S)
manifold_branch = 'L';       % Left branch of the manifold (L or R)
StableManifold = invariant_manifold(mu, Ln, manifold_ID, manifold_branch, target_orbit.Trajectory, rho, tspan);

manifold_ID = 'U';           % Unstable manifold (U or S)
manifold_branch = 'R';       % Left branch of the manifold (L or R)
UnstableManifold = invariant_manifold(mu, Ln, manifold_ID, manifold_branch, chaser_orbit.Trajectory, rho, tspan);

%% Plots
% Orbit representation
figure 
plot3(Sr(1,:),Sr(2,:),Sr(3,:),'b','LineWidth',0.9); 
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
grid on; 

figure_orbits = figure;
view(3)
hold on
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), 'b', 'LineWidth', 0.9);                         % Target's orbit
plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3), '-ob', 'LineWidth', 0.9, ...
      'MarkerIndices', floor(linspace(1,size(chaser_orbit.Trajectory,1),10)));                                                                  % Charser's initial orbit
plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth', 1.3);                                                                                               % Trasfer orbit
grid on; 
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
legend('Target orbit', 'Initial orbit', 'Transfer orbit', 'AutoUpdate', 'off')
plot3(C(1,1),C(2,1),C(3,1),'*k');                                                                                                               % Initial conditions
plot3(C(1,end),C(2,end),C(3,end),'*k');                                                                                                         % Final conditions
plot3(L(1,Ln), L(2,Ln), 0, '+k');
labels = {'$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$'};
text(L(1,Ln)-1e-3, L(2,Ln)-1e-3, 1e-2, labels{Ln});
% for i = 1:size(StableManifold.Trajectory,1)
%     ManifoldAux = shiftdim(StableManifold.Trajectory(i,:,:));
%     S = plot3(ManifoldAux(1:StableManifold.ArcLength(i),1), ManifoldAux(1:StableManifold.ArcLength(i),2), ManifoldAux(1:StableManifold.ArcLength(i),3), 'g');
%     S.Color(4) = 0.1;
% end
% for i = 1:size(UnstableManifold.Trajectory,1)
%     ManifoldAux = shiftdim(UnstableManifold.Trajectory(i,:,:));
%     U = plot3(ManifoldAux(1:UnstableManifold.ArcLength(i),1), ManifoldAux(1:UnstableManifold.ArcLength(i),2), ManifoldAux(1:UnstableManifold.ArcLength(i),3), 'r');
%     U.Color(4) = 0.1;
% end
hold off

% Propulsive acceleration plot
figure;
hold on
plot(tau, sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)*Lem/T0^2, 'k','LineWidth',1)
plot(tau, u*Lem/T0^2, 'LineWidth', 0.3)
yline(T, '--k')
xlabel('Flight time')
ylabel('$\mathbf{a}$')
legend('$a$','$a_x$','$a_y$','$a_z$')
grid on;
xlim([0 1])

figure 
hold on
plot(tau, rad2deg(atan2(u(2,:),u(1,:)))); 
hold off 
grid on;
xlabel('Time')
ylabel('$\theta$')
title('Thrust in-plane angle')

figure 
hold on
plot(tau, rad2deg(atan2(u(3,:),sqrt(u(1,:).^2+u(2,:).^2)))); 
hold off 
grid on;
xlabel('Time')
ylabel('$\phi$')
title('Thrust out-of-plane angle')
