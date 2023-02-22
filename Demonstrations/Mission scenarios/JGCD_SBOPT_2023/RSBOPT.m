%% Autonomous RVD and docking in the CR3BP  %%
% Date: 06/01/23

%% Set up
set_graphics(); 
close all

%% Trajectory generation 
% CR3BP constants 
mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
L = libration_points(mu);           % System libration points
Lem = 384400e3;                     % Mean distance from the Earth to the Moon
T0 = 28*86400;               % Characteristic time of the Earth-Moon system
Vc = 1.025e3;                       % Characteristic velocity of the Earth-Moon system

% Differential corrector set up
nodes = 10;                         % Number of nodes for the multiple shooting corrector
maxIter = 20;                       % Maximum number of iterations
tol = 1e-10;                        % Differential corrector tolerance

% Halo characteristics 
Az = 15e6;                                                          % Orbit amplitude out of the synodic plane 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 1;                                                             % Orbits around L2
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                    % Lyapunov seed parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a lyapunov orbit seed

% Correct the seed and obtain initial conditions for a lyapunov orbit
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
Az = 20e6;                                                          % Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 1;                                                             % Orbits around L1
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute a halo orbit seed 
halo_param = [1 Az Ln gamma m];                                    % Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a halo orbit seed

[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Setup of the solution method
% Numerical solver definition 
time_distribution = 'Chebyshev';        % Distribution of time intervals
basis = 'Chebyshev';                   % Polynomial basis to be use
n = [8 8 8];                        % Polynomial order in the state vector expansion
m = 200;                                % Number of sampling points
D = 2;                                 % Degree of the dynamics 

OptProblem = Problem().DefineSolver(n, basis, m, time_distribution).AddDynamics(length(n), 3, D); 

% Chaser's initial Cartesian state vector
initial_state = chaser_orbit.Trajectory(1,1:6).'; 
target_state = target_orbit.Trajectory(500,1:6).'; 

St = [target_orbit.Trajectory(500:end,1:6); target_orbit.Trajectory(1:499,1:6)];       % Target evolution

S0 = initial_state-target_state;              % Initial relative conditions
S0 = cylindrical2cartesian(S0, false);        % Initial state vector in cylindrical coordinates   
SF = zeros(6,1);                              % Final boundary conditions (rendezvous)

% Spacecraft parameters 
T = 1e-3;              % Maximum acceleration 
T = T*T0^2/Lem;          % Normalized acceleration

% Add boundary conditions
OptProblem = OptProblem.AddBoundaryConditions(S0, SF).AddParameters([mu; T; SF(2); target_orbit.Period]);

% Add functions 
OptProblem = OptProblem.AddFunctions(@(initial, final, beta, t0, tf)BoundaryConditions(initial, final, beta, t0, tf), @(params, beta, t0, tf, tau, s)ControlFunction(params, beta, t0, tf, tau, s), ...
                                     @(params, beta, t0, tf, s, u)CostFunction(params, beta, t0, tf, s, u), @(beta, P)LinConstraints(beta, P), ...
                                     @(params, beta, t0, tf, tau, s, u)NlinConstraints(params, beta, t0, tf, tau, s, u), ...
                                     @BoundsFunction, ...
                                     @(params, initial, final)InitialGuess(params, initial, final));

STM_flag = false;                   % Computation of the STM flag

%% Optimization 
iter = 1; 
time = zeros(1,iter);
options.resultsFlag = false; 
for i = 1:iter
    tic 
    [Sr, cost, u, tau] = RSBOPT_wrapper(OptProblem, St, STM_flag);
    time(i) = toc;
end

time = mean(time);

% Perfomance evaluation 
[error, merit] = figures_merit(tau, Sr(7:12,:).'); 
effort = control_effort(tau, u, false);

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
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), 'b', 'LineWidth', 0.5);           % Target's orbit
plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3), '-*b', 'LineWidth', 0.5, ...
      'MarkerIndices', floor(linspace(1,size(chaser_orbit.Trajectory,1),10)));                                                    % Charser's initial orbit
plot3(Sr(1,:)+Sr(7,:),Sr(2,:)+Sr(8,:),Sr(3,:)+Sr(9,:),'r','LineWidth', 1.5);                                                      % Trasfer orbit
grid on; 
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
legend('Target orbit', 'Initial orbit', 'Transfer orbit', 'AutoUpdate', 'off')
plot3(Sr(1,1)+Sr(7,1),Sr(2,1)+Sr(8,1),Sr(3,1)+Sr(9,1),'*k');                                                                       % Initial conditions                                                                                                             % Initial conditions
plot3(Sr(1,end)+Sr(7,end),Sr(2,end)+Sr(8,end),Sr(3,end)+Sr(9,end),'*k');                                                           % Final conditions      
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
xlabel('$t$')
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