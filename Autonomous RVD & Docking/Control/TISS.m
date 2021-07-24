%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 1: Two-impulsive rendezvous via STM %% 
% This script provides an interface to test the two-impusilve rendezvous strategy using the STM of 
% the dynamical model. 

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% The first impulse is known as targetting, aimed to nullify the relative
% position between chaser and target after some flight time tf. The second
% impulse nullifies the relative velocity between the two to complete the
% rendezvous.

% In the relative phase space, the relative particle is driven to the
% origin of the synodic frame.

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
%Time span 
dt = 1e-3;                          %Time step
tf = 0.6;                           %Rendezvous time
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
Az = 120e6;                                                 %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 5;                                                    %Number of orbits to continuate
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, period};                      %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
direction = 1;                                              %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                         %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(end,:), maxIter, tol);

%% Modelling in the synodic frame %%
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% GNC: two impulsive rendezvous, single shooting scheme %%
%Differential corrector scheme
tol = 1e-10;                                                        %Differential corrector tolerance

tic
[St, dV, state] = TISS_control(mu, tf, s0, tol, 'Position', true);  %Controller scheme
toc

%Total maneuver metrics 
effort = control_effort(tspan, dV, true);

%Error in time 
[e, merit] = figures_merit(tspan, St);

%% GNC: one impulsive rendezvous, single shooting scheme %% 
%Controller scheme
[St2, dV2, state2] = TISS_control(mu, tf, s0, tol, 'Position', false);      

%Total maneuver metrics 
effort2 = control_effort(tspan, dV2, true);

%Error in time 
[e2, merit2] = figures_merit(tspan, St2);

%% Results %% 
%Print results 
disp('SIMULATION RESULTS: ')
if (state.State)
    disp('   Two impulsive rendezvous was achieved');
    fprintf('   Delta V budget (L1 norm): %.4ei %.4ej %.4ek \n', effort(:,2));
    fprintf('   Delta V budget (L2 norm): %.4ei %.4ej %.4ek \n', effort(:,1));
else
    disp('    Two impulsive rendezvous was not achieved');
end

disp(' ');
if (state2.State)
    disp('   One impulsive rendezvous was achieved.');
    fprintf('   Delta V budget (L1 norm): %.4ei %.4ej %.4ek \n', effort2(:,2));
    fprintf('   Delta V budget (L2 norm): %.4ei %.4ej %.4ek \n', effort2(:,1));
else
     disp('   One impulsive rendezvous was not achieved.');
end

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
plot3(St(:,7), St(:,8), St(:,9)); 
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
title('Motion in the relative configuration space');

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan, St(:,7)); 
plot(tspan, St(:,8)); 
plot(tspan, St(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinates');
grid on;
legend('$x$', '$y$', '$z$');
title('Relative position in time');
subplot(1,2,2)
hold on
plot(tspan, St(:,10)); 
plot(tspan, St(:,11)); 
plot(tspan, St(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinates');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Relative velocity in time');

%Configuration space error 
figure(4)
plot(tspan, log(e)); 
xlabel('Nondimensional epoch');
ylabel('Absolute error $\log{e}$');
grid on;
title('Absolute rendezvous error in the configuration space');

figure(5)
view(3) 
hold on
c = plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'b', 'Linewidth', 0.1); 
r = plot3(St(:,7)+St(:,1), St(:,8)+St(:,2), St(:,9)+St(:,3), 'k', 'Linewidth', 0.1); 
t = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'r', 'Linewidth', 0.1);
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
text(L(1,Ln), L(2,Ln), 5e-3, '$L_2$');
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Initial orbit', 'Rendezvous arc', 'Target orbit', 'Location', 'northeast');
title('Converged rendezvous trajectory in the absolute configuration space');

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