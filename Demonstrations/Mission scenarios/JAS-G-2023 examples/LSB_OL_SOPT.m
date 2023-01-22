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
T0 = 28*86400;                      % Characteristic time of the Earth-Moon system
Vc = 1.025e3;                       % Characteristic velocity of the Earth-Moon system
n = 6;                              % Phase-space dimension

% Differential corrector set up
nodes = 10;                         % Number of nodes for the multiple shooting corrector
maxIter = 20;                       % Maximum number of iterations
tol = 1e-10;                        % Differential corrector tolerance

% Halo characteristics 
Az = 10000e3;                                                       % Orbit amplitude out of the synodic plane 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 1;                                                             % Orbits around L2
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     % Northern halo parameters
[halo_seed, ~] = object_seed(mu, halo_param, 'Halo');          % Generate a halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

% Continuate the first halo orbit to locate the chaser spacecraft
Az = 10000e3;                                                       % Orbit amplitude out of the synodic plane 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 1;                                                             % Orbits around L2
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute a halo seed 
halo_param = [-1 Az Ln gamma m];                                    % Southern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Prescribed guidance using polynomial families
% Initial conditions 
initial_state = chaser_orbit.Trajectory(1,1:6);                     % Initial chaser conditions 
target_state = target_orbit.Trajectory(800,1:6);                    % Initial target conditions
 
% Setup of the solution 
T = 0.3e-3*(T0^2/Lem);                  % Maximum thrust in non-dimensional units
GNC.Tmax = T;                           % Constrained acceleration

method = 'Minimum time'; 

% Relative solution  
iter = 25; 
Time = zeros(1,iter);
for i = 1:iter
    tic
    [Sr, u, tf, lissajous_constants] = LSB_guidance(mu, Ln, gamma, initial_state-target_state, method, 0, GNC);  
    Time(i) = toc; 
end
Time = mean(Time);

% Absolute chaser trajectory
tau = linspace(0,tf,size(Sr,2));
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
[~, Sc] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tau, target_state, options);

% Total transfer metrics 
effort = control_effort(tau, u, false)*Vc;

% Error in time 
[e, merit] = figures_merit(tau, Sr);

% Minimum and maximum control 
u_ex(1) = (Lem/T0^2)*min(sqrt(dot(u,u,1)))*1e3;
u_ex(2) = (Lem/T0^2)*max(sqrt(dot(u,u,1)))*1e3;

% Complete arc
C = Sc.'+Sr;

%% Plots
% Orbit representation
figure 
plot3(Sr(1,:), Sr(2,:), Sr(3,:), 'k', 'LineWidth', 1); 
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on; 
% yticklabels(strrep(yticklabels, '-', '$-$'));

figure_orbits = figure;
view(3)
hold on
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3),'Color','r','LineWidth', 0.9);                    % Target's orbit
plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3),'Color','b','LineWidth', 0.9);                    % Charser's initial orbit
N = plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth',1);                                                                                               % Trasfer orbit
legend('Target orbit', 'Initial orbit', 'Transfer arc', 'Location', 'northwest', 'Autoupdate', 'off');
plot3(C(1,1),C(2,1),C(3,1),'*k');                                                                                                                % Initial conditions
plot3(C(1,end),C(2,end),C(3,end),'*k');                                                                                                          % Final conditions
plot3(L(1,Ln), L(2,Ln), 0, '+k');
labels = {'$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$'};
text(L(1,Ln)-1e-3, L(2,Ln)-1e-3, 1e-2, labels{Ln});
hold off
grid on; 
% yticklabels(strrep(yticklabels, '-', '$-$'));

% Propulsive acceleration plot
figure;
hold on
plot(tau, sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)*Lem/T0^2*1e3)
xlabel('$t$')
ylabel('$\|\mathbf{u}\|$')
grid on;
% yticklabels(strrep(yticklabels, '-', '$-$'));

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
