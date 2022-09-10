%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 14/02/21
% File: test_4.m 
% Issue: 0 
% Validated: 

%% Test 4 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 4 is concerned with validating the single-parameter continuation 
% method

%% General setup 
format long; 
set_graphics();         % Set graphical environment 
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      % Integration tolerances

%% Test values and constants
% Initial conditions
mu = 0.0121505856;                                          % Reduced gravitational parameter of the system (Earth-Moon)
Lem = 384400e3;                                             % Characteristic distance of the system (Earth-Moon)
L = libration_points(mu);                                   % System libration points
Az = 45e6;                                                  % Orbit amplitude out of the synodic plane. Play with it!
Ax = 20e6;                                                  % Orbit amplitude in the synodic plane. Play with it! 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         % Normalize distances for the E-M system
Ax = dimensionalizer(Lem, 1, 1, Ax, 'Position', 0);         % Normalize distances for the E-M system
Ln = 1;                                                     % Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          % Li distance to the second primary
m = 1;                                                      % Number of periods to compute
param = [1 Az Ln gamma m];                                  % Halo orbit parameters (-1 for southern halo)

% Correction parameters 
dt = 1e-3;                                                  % Time step to integrate converged trajectories
maxIter = 20;                                               % Maximum allowed iterations in the differential correction schemes
tol = 1e-10;                                                % Differential correction tolerance 
Bif_tol = 1e-2;                                             % Bifucartion tolerance on the stability index
num = 50;                                                   % Number of orbits to continuate
direction = 1;                                              % Direction to continuate (to the Earth)
   
%% Functions
% Compute the orbits seeds 
butterfly_seed = [1.0406 0 0.1735 0 -0.0770 0];             % State vector of a butterfly orbit
axial_seed = [0.8431 0 0 0 0.1874 0.4000];                  % State vector of an axial orbit

% Halo orbit
[halo_seed, haloT] = object_seed(mu, param, 'Halo');        % Generate a halo orbit seed

% Axial orbit
[axial_orbit, ~] = differential_correction('Axis Symmetric', mu, axial_seed, maxIter, tol);

% Butterfly orbit
[butterfly_orbit, ~] = differential_correction('Plane Symmetric', mu, butterfly_seed, maxIter, tol);

% Lyapunov family 
param_lyap = [Az Az 0 0 Ln gamma m];                                  % Lyapunov orbit parameters
[lyapunov_seed, period] = object_seed(mu, param_lyap, 'Lyapunov');    % Generate a Lyapunov orbit seed

% Continuation procedure 
method = 'SPC';                                             % Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                % Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, haloT};                       % Object and characteristics to continuate
corrector = 'Plane Symmetric';                              % Differential corrector method
setup = [mu maxIter tol direction];                         % General setup

[Results_energy, ~] = continuation(num, method, algorithm, object, corrector, setup);

% Continuation procedure 
% method = 'SPC';                                             % Type of continuation method (Single-Parameter Continuation)
% algorithm = {'Period', NaN};                                % Type of SPC algorithm (on period or on energy)
% object = {'Orbit', halo_seed, haloT};                       % Object and characteristics to continuate
% corrector = ' ';                                            % Differential corrector method
% setup = [mu maxIter tol direction];                         % General setup
% 
% [Results_period, state_period] = continuation(num, method, algorithm, object, corrector, setup);

% Continuation procedure 
% method = 'PAC';                                             % Type of continuation method (Single-Parameter Continuation)
% algorithm = {'Period', NaN};                                % Type of SPC algorithm (on period or on energy)
% object = {'Orbit', halo_seed, haloT};                       % Object and characteristics to continuate
% corrector = ' ';                                            % Differential corrector method
% setup = [mu maxIter tol direction];                         % General setup
% 
% [Results_PA, state_PA] = continuation(num, method, algorithm, object, corrector, setup);

%% Plotting and results 
% Orbit plotting
figure(1) 
view(3);
hold on
for i = 1:num
  tspan = 0:dt:Results_energy.Period(i);
  [~, S] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, Results_energy.Seeds(i,:), options);
  plot3(S(:,1), S(:,2), S(:,3), 'b', 'Linewidth', 0.9);
end
plot(L(1,Ln), L(2,Ln), '+k')
text(L(1,Ln)+2e-3, L(2,Ln), '$L_1$');
hold off
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
% title('$L_{1}$ halo family');
grid on;

% %Plot results
% figure(2) 
% view(3);
% hold on
% plot3(halo_seed(:,1), halo_seed(:,2), halo_seed(:,3), 'k');
% for i = 1:num
%   tspan = 0:dt:Results_period.Period(i);
%   [~, S] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, Results_period.Seeds(i,:), options);
%   plot3(S(:,1), S(:,2), S(:,3));
% end
% hold off
% xlabel('Synodic normalized x coordinate');
% ylabel('Synodic normalized y coordinate');
% zlabel('Synodic normalized z coordinate');
% title('Converged family of orbits');
% grid on;