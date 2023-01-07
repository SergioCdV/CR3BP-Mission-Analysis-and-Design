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
T0 = 28*86400/(2*pi);               % Mean period for the Earth-Moon system
n = 6;                              % Phase-space dimension

% Differential corrector set up
nodes = 10;                         % Number of nodes for the multiple shooting corrector
maxIter = 20;                       % Maximum number of iterations
tol = 1e-10;                        % Differential corrector tolerance

% Halo characteristics 
Az = 20e6;                                                          % Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 1;                                                             % Orbits around L1
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute a halo seed 
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

butterfly_seed = [1.0406 0 0.1735 0 -0.0770 0];                     % State vector of a butterfly orbit

[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(end,:), maxIter, tol);

%% Prescribed guidance using polynomial families
% Initial conditions 
initial_state = chaser_orbit.Trajectory(1,1:n);     % Initial chaser state
target_state = target_orbit.Trajectory(500,1:n);    % Initial target state 
 
% Setup of the solution 
T = 0.5e-3*(T0^2/Lem);                  % Maximum thrust in non-dimensional units
GNC.Algorithm = 'SDRE';                 % Solver algorithm
GNC.LQR.StateMatrix = 10*eye(2);        % State error weight matrix
GNC.LQR.ControlMatrix = eye(1);         % Control effort weight matrix
GNC.Tmax = T/sqrt(3)*(T0^2/Lem);        % Constrained acceleration
GNC.TOF = 0.1*pi;                       % Maneuver time
GNC.Polynomial = 'Bernstein';           % Polynomial family to be used

method = 'Prescribed shape-based'; 

tau = linspace(0,tf,size(Sr,2));
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);

% Relative solution    
tic
for i = 1:length(tau)
    % Guidance
    [Sr, u, tf, lissajous_constants] = LSB_guidance(mu, Ln, gamma, initial_state-target_state, method, GNC.TOF, GNC);

    % Setup of the guidance scheme 
    GNC_p = 1; 

    % Propagation 
    [~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC_p), tau, [target_state initial_state-target_state], options);
end
toc 

C = Sc(:,1:6)+Sc(:,7:12);

%% Plots
% Orbit representation
figure 
plot3(Sr(1,:),Sr(2,:),Sr(3,:),'k','LineWidth',0.4); 
xlabel('Relative synodic $x$ coordinate')
ylabel('Relative synodic $y$ coordinate')
zlabel('Relative synodic $z$ coordinate')
grid on; 

figure_orbits = figure;
view(3)
hold on
xlabel('Synodic $x$ coordinate')
ylabel('Synodic $y$ coordinate')
zlabel('Synodic $z$ coordinate')
plot3(C(1,1),C(2,1),C(3,1),'*k');                                                                                                                % Initial conditions
N = plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth',0.4);                                                                                             % Trasfer orbit
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), 'LineStyle','--','Color','r','LineWidth', 0.9);  % Target's orbit
plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3),'LineStyle','-.','Color','b','LineWidth', 0.9);   % Charser's initial orbit
plot3(C(1,end),C(2,end),C(3,end),'*k');                                                                                                          % Final conditions
plot3(L(1,Ln), L(2,Ln), 0, '+k');
labels = {'$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$'};
text(L(1,Ln)-1e-3, L(2,Ln)-1e-3, 1e-2, labels{Ln});
hold off
grid on; 
legend('off')

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
