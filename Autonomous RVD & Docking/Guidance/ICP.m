%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 25/07/22 % 

%% GNC 15: Impulsive Center Manifold Phasing %% 
% This script provides an interface to test the ICP guidance core.

% The relative motion of two spacecraft in the two different halo orbit is analysed 
% and a guidance rendezvous trajectory is designed.

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
Az = 20e6;                                                          %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
chaser_orbit = target_orbit;

%% Modelling in the synodic frame %%
index = 600;                                                        %Phasing point
r_t0 = target_orbit.Trajectory(index,1:6);                          %Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);                              %Initial chaser conditions 
rho0 = r_c0-r_t0;                                                   %Initial relative conditions
s0 = [r_t0 rho0].';                                                 %Initial conditions of the target and the relative state

% Phasing 
target_orbit.Trajectory = [target_orbit.Trajectory(index:end,:); target_orbit.Trajectory(1:index,:)];

%% Generate the guidance trajectory
%Guidance trajectory
k = 3;                          % Number of phasing revolutions
dtheta = 2*pi/target_orbit.Period*(index*dt);
restriction = 'Center';
[Str, V, state(1)] = ICP_guidance(mu, Ln, gamma, target_orbit.Period, dtheta, k, [r_t0 r_c0], tol, restriction);
St = Str(:,1:6)+Str(:,7:12);
tspan = (0:size(St,1))*dt;

% Periodicity check 
target_orbit.Trajectory = repmat(target_orbit.Trajectory(1:end-1,:),k,1);
chaser_orbit.Trajectory = repmat(chaser_orbit.Trajectory(1:end-1,:),k,1);

%% Results 
%Plot results 
figure(1) 
view(3) 
hold on
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), 'b'); 
plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3), 'r'); 
plot3(St(:,1), St(:,2), St(:,3), 'g'); 
legend('Reference target orbit', 'Chaser orbit', 'Guidance orbit', 'AutoUpdate', 'off')
plot3(St(1,1), St(1,2), St(1,3), '*r'); 
plot3(target_orbit.Trajectory(1,1), target_orbit.Trajectory(1,2), target_orbit.Trajectory(1,3), '*b'); 
plot3(target_orbit.Trajectory(end,1), target_orbit.Trajectory(end,2), target_orbit.Trajectory(end,3), '*b'); 
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
title('Guidance trajectory between periodic orbits');

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan(1:size(Str,1)), Str(:,7)); 
plot(tspan(1:size(Str,1)), Str(:,8)); 
plot(tspan(1:size(Str,1)), Str(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Configuration coordinates');
grid on;
legend('$x$', '$y$', '$z$');
title('Position in time');
subplot(1,2,2)
hold on
plot(tspan(1:size(Str,1)), Str(:,10)); 
plot(tspan(1:size(Str,1)), Str(:,11)); 
plot(tspan(1:size(Str,1)), Str(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Velocity coordinates');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Velocity in time');

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
        drawnow;
        delete(T); 
        delete(V);
    end
    hold off
end
