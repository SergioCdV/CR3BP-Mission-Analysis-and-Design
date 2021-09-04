%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 25/04/21 % 

%% GNC 10: Floquet Mode Safe Control %% 
% This script provides an interface to test the Floquet mode strategy for collision avoidance maneuvers. 

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
tf = 0.6;                           %Rendezvous time
tspan = 0:dt:2*pi;                  %Integration time span

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
T = target_orbit.Period;

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                    %Number of orbits to continuate
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, T};                           %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
direction = 1;                                              %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                         %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(2,:), maxIter, tol);

%% Modelling in the synodic frame %%
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% GNC: TITA %%
%Differential corrector set up
tol = 1e-5;                                     %Differential corrector tolerance

%Cost function matrices
penalties.R = eye(3);                           %Penalty on the impulse
penalties.Q = eye(6);                           %Penalty on the state error
penalties.M  = 0.1*eye(6);                      %Penalty on the state noise

%Select measuring times 
target_points.Noise = true;                     %Boolean to account for state noise
target_points.Times = tf*rand(1,3);             %Times to measure the state noise

thruster_model.Sigma = 0.01;                    %Velocity noise dependance on the velocity impulse

%Rotational misalignment of the thrusters
thruster_model.Rotation = [1 0 0; 0 cos(pi/18) sin(pi/18); 0 -sin(pi/18) cos(pi/18)];           

%Cost function 
cost_function = 'Position';                     %Cost function to target
two_impulsive = true;                           %Two-impulsive rendezvous boolean

[St, ~, state] = TITA_control(mu, tf, s0, tol, cost_function, zeros(1,3), two_impulsive, ...
                               penalties, target_points, thruster_model);
                           
%% GNC: FMSC %% 
%Obstacle definition in space and time
index(1) = fix(0.5/dt);                         %Time location of the collision 
index(2) = fix(0.2/dt);                         %Detection time
so = [St(index(1),7:9) 0 0 0];                  %Phase space state of the object
R = 2e-3;                                       %Radius of the CA sphere
[xo, yo, zo] = sphere;                          %Collision avoidance sphere
xo = R*xo;
yo = R*yo;
zo = R*zo;

%Safety parameters 
TOC = tspan(index(1));                          %Collision time
constraint.Constrained = false;                 %No constraints on the maneuver
constraint.SafeDistance = 1e-5;                 %Safety distance at the collision time
constraint.Period = T;                          %Orbital Period

tic
[Sc, dV, tm] = FMSC_control(mu, TOC, St(index(2),1:12), 1e-5, constraint, 'Center');
toc
Sc = [St(1:index(2)-1,1:12); Sc(:,1:12)];       %Complete trajectory

%Re-insertion
if (false)
    tic
    [Ss, dVf, tm] = FMSC_control(mu, TOC, Sc(end,1:12), 1e-5, constraint, 'Stable');
    toc
    Sc = [Sc(1:end-1,1:12); Ss];                %Complete trajectory
    dV = [dV dVf];                              %Complete control law
end

ScCAM = Sc(:,1:3)+Sc(:,7:9);                    %Absolute trajectory

%Total maneuver metrics 
effort = control_effort(tspan, dV, true);
    
%% Results %% 
%Plot results 
figure(1) 
view(3) 
hold on 
surf(xo+so(1),yo+so(2),zo+so(3));
plot3(Sc(:,7), Sc(:,8), Sc(:,9), 'b'); 
plot3(St(:,7), St(:,8), St(:,9), 'r'); 
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Colliding object', 'Rendezvous and CAM arc', 'Target orbit', 'Location', 'northeast');
title('Collision avoidance trajectory in the relative configuration space');

figure(2) 
view(3) 
hold on 
surf(xo+so(1)+S(index(1),1),yo+so(2)+S(index(1),2),zo+so(3)+S(index(1),3), 'Linewidth', 0.1);
plot3(ScCAM(:,1), ScCAM(:,2), ScCAM(:,3), 'k', 'Linewidth', 0.1); 
r = plot3(St(:,1)+St(:,7), St(:,2)+St(:,8), St(:,3)+St(:,9), 'r', 'Linewidth', 0.1); 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'b', 'Linewidth', 0.1); 
plot3(S(:,1), S(:,2), S(:,3), 'r'); 
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
text(L(1,Ln)+1e-3, L(2,Ln), 5e-3, '$L_2$');
text(xo(:,1), yo(:,1), zo(:,1)+1e-3, 'Collision');
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Colliding object', 'Initial orbit', 'Rendezvous arc', 'CAM arc', 'Target orbit', 'Location', 'northeast');
title('Collision avoidance trajectory in the absolute configuration space');

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan(1:size(Sc,1)), Sc(:,7)); 
plot(tspan(1:size(Sc,1)), Sc(:,8)); 
plot(tspan(1:size(Sc,1)), Sc(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('$x$', '$y$', '$z$');
title('Relative position in time');
subplot(1,2,2)
hold on
plot(tspan(1:size(Sc,1)), Sc(:,10)); 
plot(tspan(1:size(Sc,1)), Sc(:,11)); 
plot(tspan(1:size(Sc,1)), Sc(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinates');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Relative velocity in time');