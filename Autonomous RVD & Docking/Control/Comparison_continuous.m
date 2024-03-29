%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 17/04/21 % 

%% GNC: Continuous comparison %% 
% This script provides an interface to test the differente continuous control schemes in the same mission scenario. 

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% The SMC controller is solved to drive the relative phase space vector to the origin (rendezvous condition).

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
%Phase space dimension 
n = 6; 

%Optimization 
optimization = true;                %Optimize the controller parameters 

%Time span 
dt = 1e-3;                          %Time step
tf = 2*pi;                          %Rendezvous time
tspan = 0:dt:tf;                    %Integration time span
tspann = 0:dt:2*pi;                 %Integration time span

%CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
L = libration_points(mu);           %System libration points
Lem = 384400e3;                     %Mean distance from the Earth to the Moon

%Differential corrector set up
nodes = 10;                         %Number of nodes for the multiple shooting corrector
maxIter = 20;                       %Maximum number of iterations
tol = 1e-10;                        %Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
%Halo characteristics 
Az = 50e6;                                                  %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute

%Compute a halo seed 
halo_param = [-1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Planar', mu, halo_seed, maxIter, tol);

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                    %Number of orbits to continuate
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};         %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
direction = 1;                                              %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                         %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(end,:), maxIter, tol);

%% Modelling in the synodic frame %%
index = fix(tf/dt);                                         %Rendezvous point
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                      %Reconstructed chaser motion via Encke method

%% GNC: SMC control law %%
GNC.Algorithms.Guidance = '';                   %Guidance algorithm
GNC.Algorithms.Navigation = '';                 %Navigation algorithm
GNC.Algorithms.Control = 'LQR';                 %Control algorithm
model = 'RLM';

GNC.Guidance.Dimension = 9;                     %Dimension of the guidance law
GNC.Control.Dimension = 3;                      %Dimension of the control law

GNC.System.mu = mu;                             %Systems's reduced gravitational parameter
GNC.System.Libration = [Ln gamma];              %Libration point ID

GNC.Control.LQR.Model = model;                  %LQR model
GNC.Control.SDRE.Model = model;                 %SDRE model
GNC.Control.LQR.Q = 2*eye(9);                   %Penalty on the state error
GNC.Control.LQR.M = eye(3);                     %Penalty on the control effort
GNC.Control.LQR.Reference = Sn(index,1:3);      %Penalty on the control effort
GNC.Control.SDRE.Q = 2*eye(9);                  %Penalty on the state error
GNC.Control.SDRE.M = eye(3);                    %Penalty on the control effort

%Initial conditions 
int = zeros(1,3);                               %Integral of the relative position
slqr0 = [Sn(1,:) int];                          %Initial conditions

%Compute the trajectory
tic
[~, St_lqr] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, slqr0, options);
toc 

%Error in time 
[e_lqr, merit_lqr] = figures_merit(tspan, St_lqr);

%Control law
[~, ~, u] = GNC_handler(GNC, St_lqr(:,1:6), St_lqr(:,7:end), tspan);

%Control integrals
effort_lqr = control_effort(tspan, u, false);

GNC.Algorithms.Control = 'SDRE';                %Control algorithm

%Compute the trajectory
tic
[~, St_sdre] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, slqr0, options);
toc 

%Error in time 
[e_sdre, merit_sdre] = figures_merit(tspan, St_sdre);

%Control law
[~, ~, u] = GNC_handler(GNC, St_sdre(:,1:6), St_sdre(:,7:end), tspan);

%Control integrals
effort_sdre = control_effort(tspan, u, false);

%% GNC: SMC control law %%
GNC.Algorithms.Guidance = '';               	%Guidance algorithm
GNC.Algorithms.Navigation = '';                 %Navigation algorithm
GNC.Algorithms.Control = 'SMC';                 %Control algorithm
GNC.Guidance.Dimension = 9;                     %Dimension of the guidance law
GNC.Control.Dimension = 3;                      %Dimension of the control law
GNC.System.mu = mu;                             %System reduced gravitational parameters
GNC.Control.SMC.Parameters = [1 SMC_optimization(mu, 'L1', s0, tf)];

%Re-integrate the trajectory
tic 
[~, St_smc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Error in time 
[e_smc, merit_smc] = figures_merit(tspan, St_smc);

%Control law
[~, ~, u] = GNC_handler(GNC, St_smc(:,1:6), St_smc(:,7:12), tspan);    

%Control integrals
effort_smc = control_effort(tspan, u, false);

%% Results %% 
%Plot results 
figure(1) 
view(3) 
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3)); 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
hold off
legend('Target motion', 'Chaser motion'); 
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
title('Reconstruction of the natural chaser motion');

%Plot relative phase trajectory
figure(2) 
view(3) 
plot3(St_lqr(:,7), St_lqr(:,8), St_lqr(:,9)); 
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
title('Motion in the relative configuration space');

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan, St_smc(:,7)); 
plot(tspan, St_smc(:,8)); 
plot(tspan, St_smc(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinates');
grid on;
legend('$x$', '$y$', '$z$');
title('Relative position in time');
subplot(1,2,2)
hold on
plot(tspan, St_smc(:,10)); 
plot(tspan, St_smc(:,11)); 
plot(tspan, St_smc(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinates');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Relative velocity in time');

%Configuration space error 
figure(4)
hold on
plot(tspan, log(e_lqr));
plot(tspan, log(e_sdre));
plot(tspan, log(e_smc));
hold off
xlabel('Nondimensional epoch');
ylabel('Absolute error $\log{e}$');
legend('LQR', 'SDRE', 'SMC')
grid on;
title('Absolute rendezvous error in the configuration space');

figure(5)
view(3) 
hold on
c = plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'r', 'Linewidth', 0.1); 
r = plot3(St_smc(:,7)+St_smc(:,1), St_smc(:,8)+St_smc(:,2), St_smc(:,9)+St_smc(:,3), 'b', 'Linewidth', 0.1); 
t = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'r', 'Linewidth', 0.1);
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
hold off
text(L(1,Ln)+1e-3, L(2,Ln), 0, '$L_1$');
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Target orbit', 'Rendezvous arc', 'Location', 'northeast');
title('Converged rendezvous trajectory in the absolute configuration space');

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