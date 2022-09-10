%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 25/08/22 % 

%% Center Manifold Lissajous Guidance demonstration for IAC 2022 %% 
% This script provides an interface to demonstrate the CML guidance core

% Units are non-dimensional and solutions are expressed in the synodic
% reference frame as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
% Phase space dimension 
n = 6; 

% Time span 
dt = 1e-3;                          % Time step
tf = 1;                             % Rendezvous time
tspan = 0:dt:pi;                    % Integration time span

% CR3BP constants 
mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
L = libration_points(mu);           % System libration points
Lem = 384400e3;                     % Mean distance from the Earth to the Moon

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

% Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     % Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a halo orbit seed

[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Initial conditions %%
r_t0 = target_orbit.Trajectory(100,1:6);                          % Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);                            % Initial chaser conditions 
rho0 = r_c0-r_t0;                                                 % Initial relative conditions
s0 = [r_t0 rho0].';                                               % Initial conditions of the target and the relative state

% Integration of the model
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
Sn = S;                

% Reconstructed chaser motion 
Sr = S(:,1:6)+S(:,7:12);                                          % Reconstructed chaser motion via Encke method

%% Generate the guidance trajectory
% Constraint on the trajectory
Tsyn = target_orbit.Period*chaser_orbit.Period/(target_orbit.Period+chaser_orbit.Period);
constraint.Flag = false; 
constraint.Period = Tsyn; 

% Final transfer time span
tspan = 0:dt:tf;                 

[Str, dVilg, ilg_state, S0, SM] = ILG_guidance(mu, Ln, gamma, tf, constraint, [r_t0 r_c0], tol);
Str(end,10:12) = zeros(1,3);

% Error and control valuation 
[error(:,1), merit(:,1)] = figures_merit(tspan, Str(:,7:12));
effort(:,1) = control_effort(tspan, dVilg, true);

% Comparison against the TISS controller
[Str2, dVtiss, tiss_state] = TISS_control(mu, tf, [r_t0 r_c0-r_t0].', tol, 'Position', true); 

% Error and control valuation 
[error(:,2), merit(:,2)] = figures_merit(tspan, Str2(:,7:12));
effort(:,2) = control_effort(tspan, dVtiss, true);

% Final absolute trajectories
St = Str(:,1:n);
for i = 1:size(Str,1)
    S0(i,:) = S0(:,1:n)*SM.'+S0(i,n+1:2*n);       % Natural quasi-periodic model 
    St(i,:) = Str(i,1:n)*SM.'+Str(i,n+1:2*n);     % Transfer trajectory
end

%% Results 
% Orbit plotting 
figure(1) 
view(3) 
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'b'); 
plot3(Sr(:,1), Sr(:,2), Sr(:,3), 'r'); 
plot3(St(:,1), St(:,2), St(:,3), 'k'); 
plot3(S0(:,1), S0(:,2), S0(:,3), 'g'); 
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Reference target orbit', 'Chaser orbit', 'Guidance orbit', 'Quasi-periodic guess')
title('Guidance trajectory between periodic orbits');

% Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan(1:size(Str,1)), Str(:,7)); 
plot(tspan(1:size(Str,1)), Str(:,8)); 
plot(tspan(1:size(Str,1)), Str(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Configuration coordinates');
grid on;
legend('$x$', '$y$', '$z$');
title('Position in time');
subplot(1,2,2)
hold on
plot(tspan(1:size(Str,1)), Str(:,10)); 
plot(tspan(1:size(Str,1)), Str(:,11)); 
plot(tspan(1:size(Str,1)), Str(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Velocity coordinates');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Velocity in time');

% Rendezvous animation 
if (false)
    figure(5) 
    view(3) 
    grid on;
    hold on
    plot3(Sc(1:index,1), Sc(1:index,2), Sc(1:index,3), 'k-.'); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    title('Rendezvous simulation');
    for i = 1:size(Sc,1)
        T = scatter3(Sc(i,1), Sc(i,2), Sc(i,3), 30, 'b'); 
        drawnow;
        delete(T); 
        delete(V);
    end
    hold off
end
