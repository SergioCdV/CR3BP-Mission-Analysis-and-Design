%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 17/10/22 % 

%% Continuous guidance comparison %% 
% This script provides an interface to test the differente continuous control schemes in the same mission scenario. 

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

% Optimization 
optimization = true;                % Optimize the controller parameters 

% Time span 
dt = 1e-3;                          % Time step
tf = 2*pi;                          % Rendezvous time

% CR3BP constants 
mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
L = libration_points(mu);           % System libration points
Lem = 384400e3;                     % Mean distance from the Earth to the Moon
T = 28*86400;                       % Characteritic period of the Earth-Moon system
Vc = 1.025e3;                       % Characteristic velocity of the Earth-Moon system

% Differential corrector set up
nodes = 10;                         % Number of nodes for the multiple shooting corrector
maxIter = 20;                       % Maximum number of iterations
tol = 1e-10;                        % Differential corrector tolerance

%% Halo orbit computation %%
% Halo characteristics 
Az = 10e6;                                                  % Orbit amplitude out of the synodic plane
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         % Normalize distances for the E-M system
Ln = 1;                                                     % Orbits around L1
gamma = L(end,Ln);                                          % Li distance to the second primary
m = 1;                                                      % Number of periods to compute

% Compute a halo orbit seed 
halo_param = [-1 Az Ln gamma m];                            % Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  % Generate a halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Planar', mu, halo_seed, maxIter, tol);

%% Natural motion %%
index = fix(tf/dt);                                         % Rendezvous point
r_t0 = target_orbit.Trajectory(100,1:6);                    % Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      % Initial chaser conditions 
rho0 = r_c0-r_t0;                                           % Initial relative conditions
s0 = [r_t0 rho0].';                                         % Initial conditions of the target and the relative state

% Integration of the model
tspan = 0:dt:2*pi;
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
Sn = S;                

% Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  % Reconstructed chaser motion via Encke method

%% GNC: LQR/SDRE control law %%
GNC.Algorithms.Guidance = '';                   % Guidance algorithm
GNC.Algorithms.Navigation = '';                 % Navigation algorithm
GNC.Algorithms.Control = 'LQR';                 % Control algorithm

GNC.Guidance.Dimension = 9;                     % Dimension of the guidance law
GNC.Control.Dimension = 3;                      % Dimension of the control law
GNC.Navigation.NoiseVariance = 0;               % Noise variance

GNC.System.mu = mu;                             % Systems's reduced gravitational parameter
GNC.System.Libration = [Ln gamma];              % Libration point ID

GNC.Control.LQR.Model = 'SLLM';                 % LQR model
GNC.Control.LQR.Q = 2*eye(9);                   % Penalty on the state error
GNC.Control.LQR.M = eye(3);                     % Penalty on the control effort
GNC.Control.LQR.Reference = Sn(index,1:3);      % Reference operting point

GNC.Control.SDRE.Model = 'RLM';                 % SDRE model
GNC.Control.SDRE.Q = 2*eye(9);                  % Penalty on the state error
GNC.Control.SDRE.M = eye(3);                    % Penalty on the control effort

% Augmented initial conditions 
int = zeros(1,3);                               % Integral of the relative position
s0 = [s0.' int];                                % Initial conditions

% Compute the LQR trajectory
iter = 1; 
time = zeros(1,iter);
for i = 1:iter
    tic
    [~, St_lqr] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
    time(i) = toc; 
end
Time(1) = mean(time); 

% Error in time 
[e_lqr, merit_lqr] = figures_merit(tspan, St_lqr(:,7:12));

% Control law
[~, ~, u_lqr] = GNC_handler(GNC, St_lqr(:,1:n), St_lqr(:,7:end), tspan);

% Control integrals
effort_lqr = control_effort(tspan, u_lqr, false);

GNC.Algorithms.Control = 'SDRE';                % Control algorithm

% Compute the trajectory
for i = 1:iter
    tic
    [~, St_sdre] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
    time(i) = toc; 
end
Time(2) = mean(time); 

% Error in time 
[e_sdre, merit_sdre] = figures_merit(tspan, St_sdre(:,7:12));

% Control law
[~, ~, u_sdre] = GNC_handler(GNC, St_sdre(:,1:n), St_sdre(:,7:end), tspan);

% Control integrals
effort_sdre = control_effort(tspan, u_sdre, false);

%% GNC: iLQR control law
GNC.Algorithms.Control = 'LQR';                 % Control algorithm
GNC.Control.iLQR.Mode = 'Continuous';           % iLQR solver

GNC.Control.LQR.Q = 1e3*blkdiag(eye(3), 1e-6*eye(6));           % Penalty on the state error
GNC.Control.LQR.M = eye(3);                                     % Penalty on the control effort

Tmax = 1e-3 / (4*pi^2*Lem/T^2);                                 % Maximum available acceleration
tol = [1e-4 1e-5];                                              % Convergence tolerance

% Compute the trajectory
iter = 1; 
time = zeros(1,iter);
for i = 1:iter
    tic
    [tspan_ilqr, St_ilqr, u_ilqr, ~] = iLQR_control(mu, 2, s0, GNC, Tmax, 1e-2, tol);
    time(i) = toc; 
end
Time(3) = mean(time);

% Error in time 
[e_ilqr, merit_ilqr] = figures_merit(tspan_ilqr, St_ilqr(:,7:12));

% Control integrals
effort_ilqr = control_effort(tspan_ilqr, u_ilqr, false);

%% Results %% 
% Plot natural orbit results 
figure(1) 
view(3) 
hold on
plot3(Sn(1:2000,1), Sn(1:2000,2), Sn(1:2000,3)); 
% plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
hold off
legend('Target orbit', 'Initial orbit'); 
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;

% Plot relative phase trajectory
figure(2) 
view(3) 
plot3(St_lqr(:,7), St_lqr(:,8), St_lqr(:,9)); 
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;

% Configuration space error 
figure(4)
hold on
plot(tspan, log(e_lqr));
plot(tspan, log(e_sdre));
plot(tspan_ilqr, log(e_ilqr));
hold off
xlabel('$t$');
ylabel('Absolute error $\log{e}$');
legend('LQR', 'SDRE', 'iLQR')
grid on;

figure(5)
view(3) 
hold on
annotation('arrow', [0.25 0.30], [0.75 0.78])
c = plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'k*', 'Linewidth', 0.1, 'MarkerIndices', 1:500:size(S_rc,1)); 
r = plot3(St_ilqr(:,7)+St_ilqr(:,1), St_ilqr(:,8)+St_ilqr(:,2), St_ilqr(:,9)+St_ilqr(:,3), 'b', 'Linewidth', 1); 
g = plot3(St_lqr(:,7)+St_lqr(:,1), St_lqr(:,8)+St_lqr(:,2), St_lqr(:,9)+St_lqr(:,3), 'r', 'Linewidth', 1); 
% t = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'r', 'Linewidth', 0.1);
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
hold off
text(L(1,Ln)+1e-3, L(2,Ln), 0, '$L_1$');
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;
legend('Target orbit', 'AL-iLQR', 'LQR', 'Location', 'southwest');
% axis('equal')

figure
colors = [[0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.3010 0.7450 0.9330]];
hold on
plot(tspan(:,1:2000), sqrt(dot(u_lqr(:,1:2000), u_lqr(:,1:2000),1)) * (4*pi^2*Lem/T^2) * 1e3, 'Color', colors(1,:));
plot(tspan(:,1:2000), sqrt(dot(u_lqr(:,1:2000), u_sdre(:,1:2000),1)) * (4*pi^2*Lem/T^2) * 1e3, 'Color', colors(2,:));
plot(tspan_ilqr, sqrt(dot(u_ilqr,u_ilqr,1)) * (4*pi^2*Lem/T^2) * 1e3, 'Color', colors(3,:));
hold off
grid on;
xlabel('$t$');
ylabel('$||\mathbf{u}||$');
legend('LQR', 'SDRE', 'AL-iLQR')
%%
plotTripleEvolution(tspan_ilqr, tspan, St_lqr, St_sdre, St_ilqr);

%Rendezvous animation 
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
        V = scatter3(Sc(i,1)+Sc(i,7), Sc(i,2)+Sc(i,8), Sc(i,3)+Sc(i,9), 30, 'r');

        drawnow;
        delete(T); 
        delete(V);
    end
    hold off
end

%% Auxiliary functions 
%Function to plot the three relative state evolutions in the same figure 
function plotTripleEvolution(tspan_ilqr, tspan, St_lqr, St_sdre, St_ilqr)
    %Line format array 
    lines = {'-' '-' '-.'};

    %Colors
    colors = [[0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.3010 0.7450 0.9330]];

    %Markers 
    markers = {'none', '*', 'none'};
    markers_size = [6, 5, 6];
    
    figure
    subplot(1,2,1)
    hold on
    for i = 1:2
        %Configuration space evolution
        plot(tspan(1:4000), St_lqr(1:4000,6+i), 'Color', colors(i,:), 'LineStyle', lines{1}, 'Marker', markers{1}, 'MarkerSize', markers_size(1), 'MarkerIndices', 1:800:length(tspan(1:4000)));    
    end
    legend('$x$', '$y$', 'AutoUpdate', 'off');
    for i = 1:2
        plot(tspan(1:4000), St_sdre(1:4000,6+i), 'Color', colors(i,:), 'LineStyle', lines{2}, 'Marker', markers{2}, 'MarkerSize', markers_size(2), 'MarkerIndices', 1:800:length(tspan(1:4000)));    
        plot(tspan_ilqr, St_ilqr(:,6+i), 'Color', colors(i,:), 'LineStyle', lines{3}, 'Marker', markers{3}, 'MarkerSize', markers_size(3), 'MarkerIndices', 1:20:length(tspan_ilqr)); 
    end
    hold off 
    xlabel('$t$');
    ylabel('$\mathbf{\rho}$');
    grid on;
    ax = gca;
    ax.YAxis.Exponent = 0;

    subplot(1,2,2)
    hold on
    for i = 1:2
        plot(tspan(1:4000), St_lqr(1:4000,9+i), 'Color', colors(i,:), 'LineStyle', lines{1}, 'Marker', markers{1}, 'MarkerSize', markers_size(1), 'MarkerIndices', 1:800:length(tspan(1:4000)));  
    end
    legend('$\dot{x}$', '$\dot{y}$', 'AutoUpdate', 'off');
    for i = 1:2
        plot(tspan(1:4000), St_sdre(1:4000,9+i), 'Color', colors(i,:), 'LineStyle', lines{2}, 'Marker', markers{2}, 'MarkerSize', markers_size(2), 'MarkerIndices', 1:800:length(tspan(1:4000)));    
        plot(tspan_ilqr, St_ilqr(:,9+i), 'Color', colors(i,:), 'LineStyle', lines{3}, 'Marker', markers{3}, 'MarkerSize', markers_size(3), 'MarkerIndices', 1:20:length(tspan_ilqr));        
    end
    hold off; 
    xlabel('$t$');
    ylabel('$\dot{\mathbf{\rho}}$');
    grid on;
    ax = gca;
    ax.YAxis.Exponent = 0;
end