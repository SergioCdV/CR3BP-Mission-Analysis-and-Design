%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 09/01/23
% File: LSB_families.m 
% Issue: 0 
% Validated: 

%% LSB families %%
% This scripts provides a test function to plot the studied transfer
% families in the LSB examples

%% General setup 
format long; 
set_graphics();                                             % Set graphical environment 
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      % Integration tolerances

%% Test values and constants
% Initial conditions
mu = 0.0121505856;                  % Reduced gravitational parameter of the system (Earth-Moon)
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
Ln = 2;                                                             % Orbits around L2
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute
Az = 5e6:500e3:20e6;                                                % Orbit amplitude out of the synodic plane 

figure(1)
view(3)
hold on;
for i = 1:length(Az)
    az = dimensionalizer(Lem, 1, 1, Az(i), 'Position', 0);              % Normalize distances for the E-M system
    
    % Compute a halo seed 
    halo_param = [1 az Ln gamma m];                                     % Lyapunov seed parameters
    [halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a Halo orbit seed
    
    % Correct the seed and obtain initial conditions for a lyapunov orbit
    [chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

    % Norther plot
    plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3), 'b', 'Linewidth', 0.9);
    jacobi_constant(mu, chaser_orbit.Trajectory(i,1:6).')
end
plot(L(1,Ln), L(2,Ln), '+k')
labels = {'$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$'};
text(L(1,Ln)-1e-3, L(2,Ln)-1e-3, 1e-2, labels{Ln});

%% Southern family plots
figure(1)
hold on;
for i = 1:length(Az)
    az = dimensionalizer(Lem, 1, 1, Az(i), 'Position', 0);              % Normalize distances for the E-M system
    
    % Compute a halo seed 
    halo_param = [-1 az Ln gamma m];                                    % Lyapunov seed parameters
    [halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a Halo orbit seed
    
    % Correct the seed and obtain initial conditions for a lyapunov orbit
    [chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

    % Norther plot
    plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3), 'r', 'Linewidth', 0.9);
end
hold off
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;
