%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 16/10/22 % 

%% Impulsive controllers comparison %% 
% This script provides an interface to impulsive control schemes within the same mission scenario.

% Units are non-dimensional and solutions are expressed in the synodic
% reference frame as defined by Howell, 1984.

%% Set up %%
clc 
clear 
close all; 

% Set up graphics 
set_graphics();

% Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
% Phase space dimension 
n = 6; 

% Time span 
dt = 1e-3;                          % Time step
tf = 0.6;                           % Allowed time of flight

% CR3BP constants 
mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
L = libration_points(mu);           % System libration points
Lem = 384400e3;                     % Mean distance from the Earth to the Moon
Tc = 28*86400;                      % Characteristic time of the Earth-Moon system
Vc = 1.025e3;                       % Characteristic velocity of the Earth-Moon system

% Differential corrector set up
maxIter = 100;                      % Maximum number of iterations
tol = 1e-10;                        % Differential corrector tolerance

%% Halo orbit computation %%
% Halo characteristics 
Az = 20140e3;                                                       % Orbit amplitude out of the synodic plane 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 2;                                                             % Orbits around L2
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute the halo orbit seed 
halo_param = [1 Az Ln gamma m];                                     % Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate the halo orbit seed

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
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(end,:), maxIter, tol);

%% Natural motion %%
r_t0 = target_orbit.Trajectory(100,1:6);                            % Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);                              % Initial chaser conditions 
rho0 = r_c0-r_t0;                                                   % Initial relative conditions
s0 = [r_t0 rho0].';                                                 % Initial conditions of the target and the relative state

% Integration of the model
tspan = 0:dt:2*pi; 
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
Sn = S;                

% Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                          % Reconstructed chaser motion via Encke method

%% GNC: discrete iLQR
GNC.Algorithms.Guidance = '';                   % Guidance algorithm
GNC.Algorithms.Navigation = '';                 % Navigation algorithm
GNC.Algorithms.Control = 'LQR';                 % Control algorithm

GNC.Guidance.Dimension = 9;                     % Dimension of the guidance law
GNC.Control.Dimension = 3;                      % Dimension of the control law
GNC.Navigation.NoiseVariance = 0;               % Noise variance

GNC.System.mu = mu;                             % Systems's reduced gravitational parameter
GNC.System.Libration = [Ln gamma];              % Libration point ID

GNC.Control.LQR.Model = 'SLLM';                                 % LQR model
GNC.Control.LQR.Reference = Sn(end,1:3);                        % Reference operting point
GNC.Control.LQR.Q = [eye(3) zeros(3); zeros(3) 1e-6*eye(3)];    % Penalty on the state error
% GNC.Control.LQR.Q = 2*eye(9);
GNC.Control.LQR.M = eye(3);                                     % Penalty on the control effort

% GNC.Control.SDRE.Model = 'RLM';                 % SDRE model
% GNC.Control.SDRE.Q = 2*eye(9);                  % Penalty on the state error
% GNC.Control.SDRE.M = eye(3);                    % Penalty on the control effort

GNC.Control.iLQR.Mode = 'Discrete';                                 % iLQR solver

% iLQR control law 
int = zeros(1,3);                  % Integral of the relative position
s0 = Sn(1,:);                % Initial conditions
% s0 = [Sn(1,:) int];                % Initial conditions

% Compute the homing trajectory
tol = [1e-5 1e-5];                 % Convergence tolerance

% Controller
iter = 1; 
time = zeros(1,iter);
for i = 1:iter
    tic
    [tspan_ilqr, St_ilqr, u, state{3}] = iLQR_control(mu, pi, s0, GNC, 5e-2, 5e-2, false, tol);
    [tspan_ilqr, St_ilqr_los, u, state{3}] = iLQR_control(mu, pi, s0, GNC, 5e-2, 5e-2, true, tol);    
    time(i) = toc;
end
Time(5) = mean(time);

% Error in time 
[e_ilqr, merit_ilqr] = figures_merit(tspan_ilqr, St_ilqr_los(:,7:12));

% Control integrals
effort_ilqr = control_effort(tspan_ilqr, u, true);

% Corridor angle 
p = repmat([1,1,1],size(St_ilqr,1),1);
phi(1,:) = acos(dot(p,St_ilqr_los(:,7:9),2)./(sqrt(dot(St_ilqr_los(:,7:9),St_ilqr_los(:,7:9),2)).*sqrt(dot(p,p,2))));
phi(2,:) = acos(dot(p,St_ilqr(:,7:9),2)./(sqrt(dot(St_ilqr(:,7:9),St_ilqr(:,7:9),2)).*sqrt(dot(p,p,2))));

phi = rad2deg(phi)/2;

%% Results %% 
% Plot results 
figure
view(3) 
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3)); 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
hold off
legend('Target orbit', 'Initial orbit'); 
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;

% Corridor angle 
figure
hold on
yline(25, 'k--')
plot(tspan_ilqr, phi(1,:));
plot(tspan_ilqr, phi(2,:));
hold off
xlabel('$t$');
ylabel('$\phi$');
grid on;
legend('$\beta$', 'LOS-constrained', 'unconstrained')

% Plot relative phase trajectory
figure
view(3) 
hold on
plot3(St_ilqr_los(10:end,7), St_ilqr_los(10:end,8), St_ilqr_los(10:end,9)); 
plot3(St_ilqr(10:end,7), St_ilqr(10:end,8), St_ilqr(10:end,9)); 
hold off
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
legend('LOS-constrained', 'unconstrained')
grid on;

figure
view(3) 
hold on
t = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'k', 'Linewidth', 0.8);
c = plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'k-*', 'Linewidth', 0.8, 'MarkerIndices', 1:500:size(S_rc,1));  
r = plot3(St_ilqr(:,7)+St_ilqr(:,1), St_ilqr(:,8)+St_ilqr(:,2), St_ilqr(:,9)+St_ilqr(:,3), 'Linewidth', 1); 
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
text(L(1,Ln), L(2,Ln), 5e-3, '$L_2$');
hold off
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;
legend('Target orbit', 'Initial orbit', 'AL-iLQR', 'Location', 'northwest');
 
%Rendezvous animation 
if (false)
    figure(5) 
    view(3) 
    grid on;
    hold on
    plot3(St(:,1), St(:,2), St(:,3), 'k-.'); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    title('Rendezvous simulation');
    for i = 1:size(St,1)
        T = scatter3(St(i,1), St(i,2), St(i,3), 30, 'b'); 
        V = scatter3(St(i,1)+St(i,7), St(i,2)+St(i,8), St(i,3)+St(i,9), 30, 'r');

        drawnow;
        delete(T); 
        delete(V);
    end
    hold off
end

plotTripleEvolution(tspan_ilqr, St_ilqr, St_ilqr_los)

%% Auxiliary functions 
%Function to plot the three relative state evolutions in the same figure 
function plotTripleEvolution(tspan, St_tiss, St_ilqr)
    %Line format array 
    lines = {'-' '-'};

    %Colors
    colors = [[0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.3010 0.7450 0.9330]];

    %Markers 
    markers = {'none', '*'};
    markers_size = [6, 5];
    
    figure
    subplot(1,2,1)
    hold on
    for i = 1:3
        %Configuration space evolution
        plot(tspan(30:end), St_tiss((30:end),6+i), 'Color', colors(i,:), 'LineStyle', lines{1}, 'Marker', markers{1}, 'MarkerSize', markers_size(1), 'MarkerIndices', 1:1:length(tspan((30:end)))); 
    end
    legend('$x$', '$y$', '$z$', 'AutoUpdate', 'off');
    for i = 1:3
        %Configuration space evolution
        plot(tspan((30:end)), St_ilqr((30:end),6+i), 'Color', colors(i,:), 'LineStyle', lines{2}, 'Marker', markers{2}, 'MarkerSize', markers_size(2), 'MarkerIndices', 1:10:length(tspan((30:end))));
    end
    hold off;
    xlabel('$t$');
    ylabel('$\mathbf{\rho}$');
    grid on;

    subplot(1,2,2)
    hold on
    for i = 1:3
        %Configuration space evolution
        plot(tspan((30:end)), St_tiss((30:end),9+i), 'Color', colors(i,:), 'LineStyle', lines{1}, 'Marker', markers{1}, 'MarkerSize', markers_size(1), 'MarkerIndices', 1:1:length(tspan((30:end)))); 
    end
    legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$', 'AutoUpdate', 'off');
    for i = 1:3
        %Configuration space evolution
        plot(tspan((30:end)), St_ilqr((30:end),9+i), 'Color', colors(i,:), 'LineStyle', lines{2}, 'Marker', markers{2}, 'MarkerSize', markers_size(2), 'MarkerIndices', 1:10:length(tspan((30:end))));
    end
    hold off;

    xlabel('$t$');
    ylabel('$\dot{\mathbf{\rho}}$');
    grid on;
end