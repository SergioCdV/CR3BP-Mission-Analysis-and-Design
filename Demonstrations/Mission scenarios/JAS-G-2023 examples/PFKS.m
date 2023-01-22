%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 15/01/23 % 

%% Relative Floquet Stationkeeping demonstration %% 
% This script provides an interface to test the RFSK strategies for optimal long term stationkeeping

% Units are non-dimensional and solutions are expressed in the synodic
% reference frame as defined by Howell, 1984.

%% Set up %%
% Set up graphics 
set_graphics();

% Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
% Phase space dimension 
n = 6; 

% Time span 
dt = 1e-3;                          % Time step
tf = 1.2*pi;                        % Rendezvous time

% CR3BP constants 
mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
L = libration_points(mu);           % System libration points
Lem = 384400e3;                     % Mean distance from the Earth to the Moon
T0 = 28*86400;                      % Characteristic time of the Earth-Moon system
Vc = 1.025e3;                       % Characteristic velocity of the Earth-Moon system

% Differential corrector set up
nodes = 10;                         % Number of nodes for the multiple shooting corrector
maxIter = 20;                       % Maximum number of iterations
tol = 1e-10;                        % Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
% Halo characteristics 
Az = 30e6;                                                          % Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 2;                                                             % Orbits around L1
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute a halo seed 
halo_param = [-1 Az Ln gamma m];                                    % Southern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% GNC: RFSK stationkeeping control law
% Initial conditions 
GNC.Tmax = 0.5e-3*(T0^2/Lem);

k = dimensionalizer(Lem, 1, 1, 190e3, 'Position', 0);     % Noise gain
k = [repmat(k,1,3) 4*ones(1,3)/Vc];
rho0 = normrnd(zeros(1,6),k);                             % Noisy relative initial conditions

s0 = [target_orbit.Trajectory(1,1:n) [192e3 -192e3 192e3 13 0.5 0.5]./[Lem Lem Lem Vc Vc Vc]];

% Stationkeeping trajectory
[tspan, St, u] = PFSK_def(mu, 0.01, target_orbit.Period, s0, 8, GNC);

% Error in time 
[~, merit(:,1)] = figures_merit(tspan, St(:,n+1:2*n));

% Control integrals
effort = control_effort(tspan, u, false)*Vc;
u_m(1) = Lem/T0^2*min(sqrt(dot(u,u,1)))*1E3;
u_m(2) = Lem/T0^2*max(sqrt(dot(u,u,1)))*1E3;

%% Results %% 
% Orbit plotting
figure(1) 
view(3) 
hold on
index = floor(linspace(1,size(St,1),50));
scatter3(Sn(index,1), Sn(index,2), Sn(index,3), 'b','filled'); 
plot3(St(:,1), St(:,2), St(:,3), 'k','LineWidth', 0.9); 
%plot3(Sr(:,1), Sr(:,2), Sr(:,3), 'r', 'LineWidth', 0.9); 
hold off
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;
legend('Reference orbit', 'Stationkeeping orbit')
% legend('Reference orbit', 'Stationkeeping orbit', 'Natural orbit')

% Phase-space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan, St(:,1)); 
plot(tspan, St(:,2)); 
plot(tspan, St(:,3)); 
hold off
xlabel('$t$');
ylabel('$\mathbf{\rho}$');
grid on;
legend('$x$', '$y$', '$z$');
title('Position in time');
subplot(1,2,2)
hold on
plot(tspan, St(:,4)); 
plot(tspan, St(:,5)); 
plot(tspan, St(:,6)); 
hold off
xlabel('$t$');
ylabel('$\dot{\mathbf{\rho}}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');

% Error plot
figure(4)
plot(tspan, log(e)); 
xlabel('Nondimensional epoch');
ylabel('Absolute error $\log{e}$');
grid on;
title('Absolute rendezvous error in the relative space');

figure 
plot(tspan(tspan < 0.1*target_orbit.Period), abs(alpha(tspan < 0.1*target_orbit.Period)))
xlabel('$t$')
ylabel('$\alpha_1$')
grid on; 

figure 
plot(tspan(tspan < 0.1*target_orbit.Period), log(sqrt(dot(u,u,1))*Lem/T0^2*1E3))
xlabel('$t$')
ylabel('log $||\mathbf{u}||$')
grid on;

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
