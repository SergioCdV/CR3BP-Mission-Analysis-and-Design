%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 07/02/23 % 

%% Floquet Mode Control %% 
% This script provides an interface to generate the dynamical control mission examples to demonstrate the FMSC capabilities. 

% Units are non-dimensional and solutions are expressed in the synodic reference frame as defined by Howell, 1984.

%% Set up %%
clear; 
close all; 
clc; 

% Set up graphics 
set_graphics();

% Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
% CR3BP constants 
mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
L = libration_points(mu);           % System libration points
Lem = 384400e3;                     % Mean distance from the Earth to the Moon
T0 = 28*86400/(2*pi);               % Characteristic time of the Earth-Moon system
Vc = 1.025e3;                       % Characteristic velocity of the Earth-Moon system
n = 6;                              % Phase-space dimension

% Time span 
dt = 1e-3;                          % Time step
tf = 0.6;                           % Rendezvous time
tspan = 0:dt:2*pi;                  % Integration time span

% Differential corrector set up
nodes = 10;                         % Number of nodes for the multiple shooting corrector
maxIter = 20;                       % Maximum number of iterations
tol = 1e-10;                        % Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
% Halo characteristics 
Az = 30e6;                                                  % Orbit amplitude out of the synodic plane
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         % Normalize distances for the E-M system
Ln = 1;                                                     % Orbits around L1
gamma = L(end,Ln);                                          % Li distance to the second primary
m = 1;                                                      % Number of periods to compute

% Compute a halo orbit seed 
halo_param = [-1 Az Ln gamma m];                            % Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  % Generate a halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
T = target_orbit.Period;

% Continuate the first halo orbit to locate the chaser spacecraft
Az = 30e6;                                                          % Orbit amplitude out of the synodic plane 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 1;                                                             % Orbits around L2
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute a halo seed 
halo_param = [-1 Az Ln gamma m];                                    % Southern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Natural trajectory %%
r_t0 = target_orbit.Trajectory(2,1:6);                      % Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      % Initial chaser conditions 
rho0 = r_c0-r_t0;                                           % Initial relative conditions
s0 = [r_t0 rho0].';                                         % Initial conditions of the target and the relative state

% Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
Sn = S;                

% Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  % Reconstructed chaser motion via Encke method
                 
%% GNC: FMSC %% 
% Obstacle definition in space and time
index(1) = fix(0.5/dt);                         % Time location of the collision 
index(2) = fix(0.2/dt);                         % Detection time
so = [Sn(index(1),7:9) 0 0 0];                  % Phase space state of the object

% Safety parameters 
TOC = tspan(index(1))-tspan(index(2));          % Collision time
setup.ReferencePeriod = T;                      % Orbital Period
setup.STM = 'Numerical';                        % STM model to be used
setup.Restriction = 'Center';                   % Dynamical structures to be used for the departure maneuver
setup.LibrationID = Ln;                         % Libration point of interest
setup.SafetyDistance = 100e3/Lem;               % Safety distance to the collider
setup.Reinsertion = false;                      % Do not re-insert into orbit
tic
[tspan, Sc, dV, state(1)] = FMSC_control(mu, [TOC 2 * TOC], Sn(index(2),1:12), setup);
toc

% Distance to the object
d(1) = norm(Sc(end,7:9))*Lem;

tol = [1e-7 1e-3];            % Convergence toleranes 
N = 50;                       % Number of impulses
method = 'Primal';            % Solver method    
integrator = 'Numerical';     % Integrator to be used

% Controller scheme
tic
[tspan_misg, St_misg, dV_misg, state(2)] = MISG_control(mu, Ln, 3 * TOC, Sc(end,1:2*n).', method, integrator, N, tol);  
toc;

dV = [dV(:,1:end-1) dV_misg];
tspan = [tspan(1:end-1) tspan(end)+tspan_misg];
Sc = [Sc(1:end-1,:); St_misg];

% Total maneuver metrics 
effort = Vc*control_effort(tspan, dV, true);

% Forward the trajectory
Sc = [Sn(1:index(2)-1,1:12); Sc(:,1:12)];       % Complete trajectory
ScCAM = Sc(:,1:6)+Sc(:,7:12);                   % Absolute trajectory

tspan = [t(1:index(2)-1).' t(index(2)-1).'+tspan];
    
%% Results %% 
% Plot results 
figure(1) 
view(3) 
hold on 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'b'); 
plot3(S(:,1), S(:,2), S(:,3), 'r'); 
r = plot3(Sn(:,1)+Sn(:,7), Sn(:,2)+Sn(:,8), Sn(:,3)+Sn(:,9), 'r');
plot3(ScCAM(:,1), ScCAM(:,2), ScCAM(:,3), 'k'); 
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
scatter3(so(1)+S(index(1),1), so(2)+S(index(1),2), so(3)+S(index(1),3), 'k', 'filled');
text(L(1,Ln)+1e-3, L(2,Ln), 5e-3, '$L_2$');
hold off
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;
legend('Initial orbit', 'Target orbit', 'Rendezvous arc', 'CAM arc');
% yticklabels(strrep(yticklabels, '-', '$-$'));

% Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan(1:size(Sc,1)), Sc(:,7)); 
plot(tspan(1:size(Sc,1)), Sc(:,8)); 
plot(tspan(1:size(Sc,1)), Sc(:,9)); 
hold off
xlabel('$t$');
ylabel('$\rho$');
grid on;
legend('$x$', '$y$', '$z$');
subplot(1,2,2)
hold on
plot(tspan(1:size(Sc,1)), Sc(:,10)); 
plot(tspan(1:size(Sc,1)), Sc(:,11)); 
plot(tspan(1:size(Sc,1)), Sc(:,12)); 
hold off
xlabel('$t$');
ylabel('$\dot{\rho}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
% yticklabels(strrep(yticklabels, '-', '$-$'));

% Rendezvous animation 
if (false)
    figure(6) 
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