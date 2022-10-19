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

%% GNC: two impulsive rendezvous, single shooting scheme %%
% Differential corrector scheme setup
tspan = 0:dt:tf;                        % Manoeuvre integration time span
tol = 1e-10;                            % Differential corrector tolerance

% Controller scheme
iter = 1; 
time = zeros(1,iter);
for i = 1:iter
    tic
    [St_tiss, dV, state{1}] = TISS_control(mu, tf, s0, tol, 'Position', true);  
    time(i) = toc;
end
Time(1) = mean(time);

% Total maneuver metrics 
effort_tiss = control_effort(tspan, dV, true);

% Error in time 
[e_tiss, merit_tiss] = figures_merit(tspan, St_tiss(:,7:12));

%% GNC: multi impulsive rendezvous, single shooting scheme %%
% Differential corrector set up
tol = 1e-10;                                  % Differential corrector tolerance

% Select impulsive times 
%times = tf*rand(1,6);                        % Times to impulse the spacecraft
times = [0.5743 0.2912 0.4802 0.0851 0.2531 0.5494];

% Compute the control law
impulses.Number = length(times);              % Number of impulses
impulses.Weights = eye(impulses.Number*3);    % Weightening matrix
impulses.Times = times;                       % Impulses times
cost = 'Position';                            % Cost function to target

% Controller scheme
time = zeros(1,iter);
for i = 1:iter
    tic
    [St_miss, dV, state{2}] = MISS_control(mu, tf, s0, tol, cost, impulses);
    time(i) = toc;
end
Time(2) = mean(time);

% Control effort 
effort_miss = control_effort(tspan, dV, true);

% Error in time 
[e_miss, merit_miss] = figures_merit(tspan, St_miss(:,7:12));

%% GNC: multi-impulsive staging rendezvous
% Scheme setup
tol = [1e-7 1e-3];             % Convergence toleranes 
N = 60;                        % Number of impulses
method = 'MPC';                % Solver method                                

% Controller scheme
time = zeros(1,iter);
for i = 1:iter
    tic
    [tspan_misg, St_misg, dV, state{3}] = MISG_control(mu, tf, s0, method, N, tol);  
    time(i) = toc;
end
Time(3) = mean(time);

% Total maneuver metrics 
effort_misg = control_effort(tspan_misg, dV, true);

% Error in time 
[e_misg, merit_misg] = figures_merit(tspan_misg, St_misg(:,7:12));

%% GNC: MPC multi impulsive rendezvous %%
% Set up of the optimization
method = 'NPL';                               % Method to solve the problem
core = 'Linear';                              % Number of impulses
TOF = tf;                                     % Time of flight
cost_function = 'Position';                   % Target a position rendezvous

% Thruster characteristics 
Tmin = -1e-3;                                 % Minimum thrust capability (in velocity impulse)
Tmax = 1e-3;                                  % Maximum thrust capability (in velocity impulse)

% Main computation 
iter = 1; 
time = zeros(1,iter);
for i = 1:1
    tic
    [tspan_mpc, St_mpc, dV, state{4}] = MPC_control(mu, cost_function, Tmin, Tmax, TOF, s0, core, method);
    time(i) = toc;
end
Time(4) = mean(time);
%%
tspan_mpc = 0:1e-2:TOF;
[~, Sopt] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, [s0 reshape(eye(6), [1 6^2])], options);
iter = 1; 
for i = 1:iter
    tic
    [dV] = OPTI_guidance(cost_function, Tmin, Tmax, tspan, Sopt, 'Corrector', cost_function);
    time(i) = toc;
end
Time(4) = mean(time);

for i = 1:length(tspan_mpc)-1
    s0opt = Sopt(i,:); 
    s0opt(10:12) = s0opt(10:12)+dV(:,i).';
   [~, aux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), [0 tspan_mpc(2)-tspan_mpc(1)], s0opt, options);
   Sopt(i+1,:) = aux(end,:);
end

% Control integrals
effort_mpc = control_effort(tspan_mpc, dV, true);

% Error in time 
[e_mpc, merit_mpc] = figures_merit(tspan_mpc, St_mpc(:,7:12));

%% GNC: discrete iLQR
GNC.Control.iLQR.Mode = 'Discrete';                                 % iLQR solver

GNC.Control.LQR.Q = [eye(3) zeros(3); zeros(3) 1e-6*eye(3)];    % Penalty on the state error
% GNC.Control.LQR.Q = blkdiag(eye(3), 1e-4*eye(6));
GNC.Control.LQR.M = eye(3);                                         % Penalty on the control effort

% iLQR control law 
int = zeros(1,3);                  % Integral of the relative position
s0 = [Sn(1,:)];                    % Initial conditions

% Compute the trajectory
tol = [1e-4 1e-5];                 % Convergence tolerance

% Controller
time = zeros(1,iter);
for i = 1:iter
    tic
    [tspan_ilqr, St_ilqr, u, state{3}] = iLQR_control(mu, pi, s0, GNC, 5e-2, 5e-2, tol);    
    time(i) = toc;
end
Time(5) = mean(time);
 
% Error in time 
[e_ilqr, merit_ilqr] = figures_merit(tspan_ilqr, St_ilqr(:,7:12));

% Control integrals
effort_ilqr = control_effort(tspan_ilqr, u, true);

%% Results %% 
% Plot results 
figure(1) 
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

% Plot relative phase trajectory
figure(2) 
view(3) 
plot3(St_miss(:,7), St_miss(:,8), St_miss(:,9)); 
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;

% Configuration space error 
figure(3)
hold on 
plot(tspan, log(e_tiss));
plot(tspan, log(e_miss));
plot(tspan, log(e_opt)); 
plot(tspan_misg, log(e_misg)); 
plot(tspan_ilqr, log(e_ilqr)); 
xlabel('$t$');
ylabel('Absolute error $\log{e}$');
legend('TI', 'MI', 'Opt-MI', 'MISG', 'AL-iLQR'); 
grid on;
axes('position', [.22 .47 .6 .25])
box on
indexOfInterest = (tspan > 0.1) & (tspan < 0.5); 
hold on
plot(tspan(indexOfInterest), e_tiss(indexOfInterest))  
plot(tspan(indexOfInterest), e_miss(indexOfInterest))  
plot(tspan(indexOfInterest), e_opt(indexOfInterest))  
plot(tspan_misg(indexOfInterest), e_misg(indexOfInterest)) 
plot(tspan_ilqr(indexOfInterest), e_ilqr(indexOfInterest)) 
hold off
axis tight

figure(4)
view(3) 
hold on
c = plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'b', 'Linewidth', 0.8); 
r = plot3(St_miss(:,7)+St_miss(:,1), St_miss(:,8)+St_miss(:,2), St_miss(:,9)+St_miss(:,3), 'k', 'Linewidth', 0.8); 
t = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'r', 'Linewidth', 0.8);
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
text(L(1,Ln), L(2,Ln), 5e-3, '$L_2$');
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Initial orbit', 'Rendezvous arc', 'Target orbit', 'Location', 'northwest');
 
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

plotTripleEvolution(tspan, St_mpc, St_tiss, St_miss);

%% Auxiliary functions 
%Function to plot the three relative state evolutions in the same figure 
function plotTripleEvolution(tspan, St_mpc, St_tiss, St_miss)
    %Assemble the three of them in a single array 
    St = [St_mpc; St_tiss; St_miss]; 

    %Line format array 
    lines = {'-' '-.' '-'};

    %Colors
    colors = [[0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.3010 0.7450 0.9330]];

    %Markers 
    markers = {'none', 'none', '*'};
    markers_size = [6, 5, 6];
    
    figure
    for i = 1:3
        %Configuration space evolution
        subplot(1,2,1)
        hold on
        plot(tspan, St((i-1)*length(tspan)+1:i*length(tspan),7), 'Color', colors(1,:), 'LineStyle', lines{i}, 'Marker', markers{i}, 'MarkerSize', markers_size(i), 'MarkerIndices', 1:80:length(tspan)); 
        plot(tspan, St((i-1)*length(tspan)+1:i*length(tspan),8), 'Color', colors(2,:), 'LineStyle', lines{i}, 'Marker', markers{i}, 'MarkerSize', markers_size(i), 'MarkerIndices', 1:80:length(tspan)); 
        plot(tspan, St((i-1)*length(tspan)+1:i*length(tspan),9), 'Color', colors(3,:), 'LineStyle', lines{i}, 'Marker', markers{i}, 'MarkerSize', markers_size(i), 'MarkerIndices', 1:80:length(tspan)); 
        hold off
    end
    legend('$x$', '$y$', '$z$');
    xlabel('Nondimensional epoch');
    ylabel('Relative configuration coordinates');
    grid on;
    title('Relative position in time');

    for i = 1:3
        subplot(1,2,2)
        hold on
        plot(tspan, St((i-1)*length(tspan)+1:i*length(tspan),10), 'Color', colors(1,:), 'LineStyle', lines{i}, 'Marker', markers{i}, 'MarkerSize', markers_size(i), 'MarkerIndices', 1:80:length(tspan)); 
        plot(tspan, St((i-1)*length(tspan)+1:i*length(tspan),11), 'Color', colors(2,:), 'LineStyle', lines{i}, 'Marker', markers{i}, 'MarkerSize', markers_size(i), 'MarkerIndices', 1:80:length(tspan)); 
        plot(tspan, St((i-1)*length(tspan)+1:i*length(tspan),12), 'Color', colors(3,:), 'LineStyle', lines{i}, 'Marker', markers{i}, 'MarkerSize', markers_size(i), 'MarkerIndices', 1:80:length(tspan)); 
        hold off
        legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$', 'AutoUpdate', 'off');
    end
    xlabel('Nondimensional epoch');
    ylabel('Relative velocity coordinates');
    grid on;
    title('Relative velocity in time');
end