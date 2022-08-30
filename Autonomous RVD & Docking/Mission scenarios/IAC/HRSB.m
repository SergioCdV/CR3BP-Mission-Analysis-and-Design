%% Autonomous RVD and docking in the CR3BP  %%
% Date: 27/08/22

%% Halo raising demonstration for IAC 2022 %% 
% This script provides an interface to demonstrate the HRSB guidance core.

% The relative motion of two spacecraft in the same halo orbit around L1 in the
% Earth-Moon system is analyzed.

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.


%% Set up
set_graphics(); 
close all

%% Trajectory generation 
%CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
L = libration_points(mu);           %System libration points
Lem = 384400e3;                     %Mean distance from the Earth to the Moon
T0 = 28*86400/(2*pi);               %Mean period for the Earth-Moon system

%Differential corrector set up
nodes = 10;                         %Number of nodes for the multiple shooting corrector
maxIter = 20;                       %Maximum number of iterations
tol = 1e-10;                        %Differential corrector tolerance

%Halo characteristics 
Az = 20e6;                                                          %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                                     %Bifucartion tolerance on the stability index
num = 5;                                                            %Number of orbits to continuate
method = 'SPC';                                                     %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                        %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};                 %Object and characteristics to continuate
corrector = 'Plane Symmetric';                                      %Differential corrector method
direction = 1;                                                      %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                                 %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);

%Halo characteristics 
Az = 15e6;                                                          %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Setup of the solution method
time_distribution = 'Chebyshev';        % Distribution of time intervals
basis = 'Chebyshev';                    % Polynomial basis to be use
dynamics = 'Euler';                     % Dynamics parametrization to be used
n = [10 10 10];                         % Polynomial order in the state vector expansion
m = 500;                                % Number of sampling points
cost_function = 'Minimum energy';       % Cost function to be minimized

% System data 
system.mu = mu;     
system.Time = T0; 
system.Distance = Lem; 

% Spacecraft propulsion parameters 
T = 5e-2;     % Maximum acceleration 
K = 0;        % Initial input revolutions 

% Setup 
options.STM = false; 
options.order = n; 
options.basis = basis;
options.grid = time_distribution; 
options.nodes = m; 
options.formulation = dynamics; 
options.cost_function = cost_function;
options.resultsFlag = false; 
options.animations = false;  

%% Results
% Chaser state evolution setup
dt = 1e-3;                                              % Time step
tspan = (0:dt:chaser_orbit.Period).';                   % Integration time for the original chaser orbit
chaser.Trajectory = [tspan chaser_orbit.Trajectory];    % Chaser evolution

% Compute the relative orbit 
Tr = max(chaser_orbit.Period, target_orbit.Period);
rho0 = [chaser_orbit.Trajectory(1,1:6) target_orbit.Trajectory(1e3,1:6)-chaser_orbit.Trajectory(1,1:6)];
tspan = 0:dt:Tr;
[~, Sr] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, rho0, options);

% Simple solution    
theta = linspace(0,2*pi,2);
dV = zeros(size(theta));

options.nodes = 500; 
for i = 1:length(theta)
    tic
    [~, dV(i), ~, ~, ~, ~, ~, ~] = hrsb_optimization(system, target_orbit.Trajectory(:,1:6), chaser, theta(i), K, T, options);
    toc 
end

options.nodes = 500; 
options.resultsFlag = true;
[~, index] = sort(dV); 
tic
[C, dV, u, tf, tfapp, tau, exitflag, output] = hrsb_optimization(system, target_orbit.Trajectory(:,1:6), chaser, theta(index(1)), K, T, options);
toc 

% Average results 
iter = 0; 
time = zeros(1,iter);
options.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, tf, tfapp, tau, exitflag, output] = hrsb_optimization(system, target_orbit.Trajectory, chaser, K, T, options);
    time(i) = toc;
end

time = mean(time);

if (options.STM)
    % Analysis of the STM 
    d = zeros(1,size(C,2));         % STM determinant
    alpha = zeros(6,size(C,2));
    
    for i = 1:size(C,2)
        STM = reshape(C(end-35:end,i),[6 6]); 
        d(i) = det(STM);
    
        [V, lambda] = eig(STM); 
        alpha(:,i) = V^(-1)*C(1:6,i);
    end
end

%% Manifolds computation
rho = 1;                     % Number of manifold fibers to compute
tspan = 0:dt:0.1*tf;         % Integration timespan

manifold_ID = 'S';           % Stable manifold (U or S)
manifold_branch = 'L';       % Left branch of the manifold (L or R)
StableManifold = invariant_manifold(mu, Ln, manifold_ID, manifold_branch, target_orbit.Trajectory, rho, tspan);

manifold_ID = 'U';           % Unstable manifold (U or S)
manifold_branch = 'R';       % Left branch of the manifold (L or R)
UnstableManifold = invariant_manifold(mu, Ln, manifold_ID, manifold_branch, chaser_orbit.Trajectory, rho, tspan);

%% Plots
figure_orbits = figure;
view(3)
hold on
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), '--b', 'LineWidth', 0.4);                         % Target's orbit
plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3), '--ob', 'LineWidth', 0.4, ...
      'MarkerIndices', floor(linspace(1,size(chaser_orbit.Trajectory,1),10)));                                                                  % Charser's initial orbit
plot3(C(1,:),C(2,:),C(3,:),'r','LineWidth', 1);                                                                                                 % Trasfer orbit
grid on; 
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
legend('Reference target orbit', 'Chaser orbit', 'Guidance transfer orbit', 'AutoUpdate', 'off')
plot3(C(1,1),C(2,1),C(3,1),'*r');                                                                                                               % Initial conditions
plot3(C(1,end),C(2,end),C(3,end),'*r');                                                                                                         % Final conditions
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