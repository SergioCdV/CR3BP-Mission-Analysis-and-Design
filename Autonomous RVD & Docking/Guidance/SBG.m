%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 14/07/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif

%% Trajectory generation 
%CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
L = libration_points(mu);           %System libration points
Lem = 384400e3;                     %Mean distance from the Earth to the Moon

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
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(end,:), maxIter, tol);

% %Halo characteristics 
% Az = 20e6;                                                          %Orbit amplitude out of the synodic plane. 
% Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
% Ln = 1;                                                             %Orbits around L1
% gamma = L(end,Ln);                                                  %Li distance to the second primary
% m = 1;                                                              %Number of periods to compute
% 
% %Compute a halo seed 
% halo_param = [-1 Az 2 L(end,2) m];                                     %Northern halo parameters
% [halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed
% 
% %Correct the seed and obtain initial conditions for a halo orbit
% [chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Setup of the solution method
time_distribution = 'Linear';           % Distribution of time intervals
basis = 'Chebyshev';                    % Polynomial basis to be use
n = [10 10 10];                         % Polynomial order in the state vector expansion
m = 100;                                % Number of sampling points

mu = 0.0121505;                         % Earth-Moon reduced gravitational parameter
L = libration_points(mu);               % System libration points
Lem = 384400e3;                         % Mean distance from the Earth to the Moon
T0 = 28*86400/(2*pi);                   % Mean period for the Earth-Moon system

system.mu = mu;     
system.Time = T0; 
system.Distance = Lem; 

% Chaser's initial Cartesian state vector
initial_state = chaser_orbit.Trajectory(1,1:6); 

% Target's final Cartesian state vector
final_state = target_orbit.Trajectory(1000,1:6); 

% Spacecraft propulsion parameters 
T = 0.05e-3;     % Maximum acceleration 

% Initial input revolutions 
K = 0;

% Setup 
options.resultsFlag = true; 
options.animations = false; 

%% Results
% Simple solution    
tic
[C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_state, final_state, K, T, m, time_distribution, basis, n, options);
toc 
S = cylindrical2cartesian(C,true);

% Average results 
iter = 0; 
time = zeros(1,iter);
options.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, time_distribution, basis, n, options);
    time(i) = toc;
end

time = mean(time);

%% Plots
% Orbit representation
figure_orbits = figure;
view(3)
hold on
xlabel('Synodic $x$ coordinate')
ylabel('Synodic $y$ coordinate')
zlabel('Synodic $z$ coordinate')
plot3(S(1,1),S(2,1),S(3,1),'*k');                                                                                                                % Initial conditions
plot3(S(1,:),S(2,:),S(3,:),'k','LineWidth',1);                                                                                                   % Chaser's orbit
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), 'LineStyle','--','Color','r','LineWidth',0.3);   % Target's orbit
plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3),'LineStyle','-.','Color','b','LineWidth',0.3);    % Charser's initial orbit
plot3(S(1,end),S(2,end),S(3,end),'*k');                                                                                                          % Final conditions
plot3(L(1,Ln), L(2,Ln), 0, '+k');
labels = {'$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$'};
text(L(1,Ln)-1e-3, L(2,Ln)-1e-3, 1e-2, labels{Ln});
hold off
grid on; 
legend('off')

% Propulsive acceleration plot
figure;
%title('Spacecraft acceleration in time')
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
