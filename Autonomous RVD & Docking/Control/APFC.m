%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 12: APF guidance-control law %% 
% This script provides an interface to test Artificial Potential Functions guidance laws for
% rendezvous missions.

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

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

%Time span 
dt = 1e-3;                          %Time step
tf = pi;                            %Rendezvous time
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
Az = 60e6;                                                          %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 2;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                                     %Bifucartion tolerance on the stability index
num = 2;                                                            %Number of orbits to continuate
method = 'SPC';                                                     %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                        %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};                 %Object and characteristics to continuate
corrector = 'Plane Symmetric';                                      %Differential corrector method
direction = 1;                                                      %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                                 %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(2,:), maxIter, tol);

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

%% APF guidance scheme
%Compute some random objects  in the relative phase space 
So(:,1) = [-0.0078 -0.0097 0.0167].';
So(:,2) = [0.001210073304029 0.000636034414309 -0.000232007357106].';

%Controller scheme penalties
Penalties.AttractivePenalty = eye(3);            %Penalty on the distance to the origin
Penalties.RepulsivePenalty = eye(3);             %Penalty on the distance to the obstacles 
Penalties.RepulsiveWidth = 1e-3;                 %Repulsive width
Penalties.Gain = 2;                              %Maneuver gain                        

%Safety message 
safe_corridor.Safety = true;                     %Apply safe corridor constraints
safe_corridor.Parameters(1) = deg2rad(10);       %Safety corridor angle
safe_corridor.Parameters(2) = 0;                 %Safety distance to the docking port
safe_corridor.Parameters(3:4) = [1.5 1];         %Dimensions of the safety corridor

%Compute the guidance law
tic
[St, u, state] = APF_control(mu, safe_corridor, Penalties, So, tf, s0);
toc 

%Performance indices
effort = control_effort(tspan, u, true);         %Control effort made
[e, merit] = figures_merit(tspan, St);           %Figures of merit of the rendezvous error 

%% Results 
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
hold on
for i = 1:size(So,2)
    scatter3(So(1,i), So(2,i), So(3,i), 'k', 'filled');
end
plot3(St(:,7), St(:,8), St(:,9)); 
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
title('Motion in the relative configuration space');

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan(1:2000), St(1:2000,7)); 
plot(tspan(1:2000), St(1:2000,8)); 
plot(tspan(1:2000), St(1:2000,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinates');
grid on;
legend('$x$', '$y$', '$z$');
title('Relative position in time');
subplot(1,2,2)
hold on
plot(tspan(1:2000), St(1:2000,10)); 
plot(tspan(1:2000), St(1:2000,11)); 
plot(tspan(1:2000), St(1:2000,12)); 
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
c = plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'b'); 
r = plot3(St(:,7)+St(:,1), St(:,8)+St(:,2), St(:,9)+St(:,3), 'k'); 
t = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'r');
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
hold off
text(L(1,Ln)+1e-3, L(2,Ln), 5e-20, '$L_2$');
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Initial orbit', 'Rendezvous arc', 'Target orbit');
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