%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 19/12/21 % 

%% ISTS Conference paper: Dynamical Control %% 
% This script provides an interface to generate the dynamical control mission examples for the 31st ISTS. 

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

%% GNC: RTPS %%
%Differential corrector set up
tol = 1e-8;                                   %Differential corrector tolerance

%Select impulsive times 
times = [0 tf*rand(1,5)];                     %Times to impulse the spacecraft

%Compute the control law
impulses.Number = length(times);              %Number of impulses
impulses.Weights = eye(impulses.Number*3);    %Weightening matrix
impulses.Times = times;                       %Impulses times

cost = 'Position';                            %Cost function to target

%Controller scheme
tic
[St, ~, state] = MISS_control(mu, tf, s0, tol, cost, impulses);
toc 
                           
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
TOC = tspan(index(1))-tspan(index(2));          %Collision time
constraint.Constrained = false;                 %No constraints on the maneuver
constraint.SafeDistance = 1e-5;                 %Safety distance at the collision time
constraint.Period = T;                          %Orbital Period
constraint.Energy = true;                       %Energy constraint

tic
[Sc, dV, tm] = FMSC_control(mu, TOC, St(index(2),1:12), 1e-5, constraint, 'Center');
toc
[~, Sc2] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), 0:dt:0.3, Sc(end,1:12), options);
Sc = [Sc; Sc2];
Sc = [St(1:index(2)-1,1:12); Sc(:,1:12)];       %Complete trajectory
ScCAM = Sc(:,1:6)+Sc(:,7:12);                   %Absolute trajectory

%Total maneuver metrics 
effort = control_effort(tspan(1:size(Sc,1)), dV, true);

%% GNC: MLQR control law
Sg = jacobi_constant(mu, St(end,1:n).');            %Reference Jacobi Constant

GNC.Algorithms.Guidance = '';                       %Guidance algorithm
GNC.Algorithms.Navigation = '';                     %Navigation algorithm
GNC.Algorithms.Control = 'MLQR';                    %Control algorithm

GNC.Guidance.Dimension = 9;                         %Dimension of the guidance law
GNC.Control.Dimension = 3;                          %Dimension of the control law

GNC.System.mu = mu;                                 %Systems's reduced gravitational parameter
GNC.Control.MLQR.Q = eye(2);                        %Penalty on the state error
GNC.Control.MLQR.M = eye(3);                        %Penalty on the control effort
GNC.Control.MLQR.Reference = Sg;                    %Penalty on the control effort
GNC.Control.MLQR.Period = target_orbit.Period;      %Penalty on the control effort

%Floquet mode matrix
s0 = [target_orbit.Trajectory(1,1:n) reshape(eye(6), [1 n^2])];
[~, Snaux] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), 0:dt:target_orbit.Period, s0, options);

[E, lambda] = eig(reshape(Snaux(end,n+1:end), [n n]));
GNC.Control.MLQR.FloquetModes = diag(log(diag(lambda))/target_orbit.Period);
for i = 1:size(E,2)
    E(:,i) = E(:,i)/lambda(i,i);
end
GNC.Control.MLQR.FloquetDirections = E; 

%Initial conditions 
k = 1e-7;                                       %Noise gain
r_t0 = St(end,1:n);                             %Initial guidance target conditions
s0 = r_t0+k*rand(1,6);                          %Noisy initial conditions
tspann = 0:dt:2*pi;                             %New integration time span

%Compute the trajectory
s0 = [r_t0 s0-r_t0 reshape(eye(n), [1 n^2])];
[~, Sr] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspann, s0, options);
tic
[~, Staux1] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspann(tspann < 0.1), s0, options);
[~, Staux2] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspann(size(Staux1,1):end), Staux1(end,:), options);
toc
Stsk = [Staux1; Staux2(2:end,:)];

%Error in time 
[e, merit] = figures_merit(tspann, [Stsk(:,1:n) abs(Stsk(:,1:n)-Sn(1:size(Stsk,1),1:n))]);

%Control law
[~, ~, u] = GNCc_handler(GNC, Staux1(:,1:n), Staux1(:,n+1:end), tspann(tspann < 0.1));

%Control integrals
energy = control_effort(tspann(tspann < 0.1), u, false);

%Final absolute trajectory
Stsk = Stsk(:,1:n)+Stsk(:,n+1:2*n);
Sr = Sr(:,1:n)+Sr(:,n+1:2*n);
    
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
legend('Colliding object', 'Rendezvous and CAM arc', 'Target orbit');
title('Collision avoidance trajectory in the relative configuration space');

figure(2) 
view(3) 
hold on 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'b'); 
plot3(S(:,1), S(:,2), S(:,3), 'r'); 
r = plot3(St(:,1)+St(:,7), St(:,2)+St(:,8), St(:,3)+St(:,9), 'r', 'Linewidth', 0.1);
plot3(ScCAM(:,1), ScCAM(:,2), ScCAM(:,3), 'k'); 
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
scatter3(so(1)+S(index(1),1), so(2)+S(index(1),2), so(3)+S(index(1),3), 'k', 'filled');
text(L(1,Ln)+1e-3, L(2,Ln), 5e-3, '$L_2$');
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Initial orbit', 'Target orbit', 'Rendezvous arc', 'CAM arc');
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

%Stationkeeping analysis 
figure(5) 
view(3) 
hold on 
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), 'b'); 
plot3(Stsk(:,1), Stsk(:,2), Stsk(:,3), 'k'); 
plot3(Sr(:,1), Sr(:,2), Sr(:,3), 'r'); 
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
text(L(1,Ln)+1e-3, L(2,Ln), 5e-3, '$L_2$');
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Target orbit', 'Stationkeeping orbit', 'Natural orbit', 'Location', 'northeast');
title('Stationkeeping trajectory in the absolute configuration space');

%Rendezvous animation 
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