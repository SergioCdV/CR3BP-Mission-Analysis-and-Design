%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 25/08/22 % 

%% Center Manifold Quasi-Periodic Phasing %% 
% This script provides an interface to demonstrate the QPP guidance core.

% The relative motion of two spacecraft in the two different halo orbit is analysed 
% and a guidance rendezvous trajectory is designed.

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Set up %%
clear; 
close all; 

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
Az = 60e6;                                                          %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Planar', mu, halo_seed, maxIter, tol);
chaser_orbit = target_orbit;

%% Modelling in the synodic frame %%
index = 5e2;                                                        %Phasing point
r_t0 = target_orbit.Trajectory(index,1:6);                          %Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);                              %Initial chaser conditions 
rho0 = r_c0-r_t0;                                                   %Initial relative conditions
s0 = [r_t0 rho0].';                                                 %Initial conditions of the target and the relative state

% Phasing 
target_orbit.Trajectory = [target_orbit.Trajectory(index:end,:); target_orbit.Trajectory(1:index,:)];

%% Generate the guidance trajectory
% Guidance trajectory parameters
k = 10;                                         % Number of phasing revolutions
dtheta = 2*pi/target_orbit.Period*(index*dt);   % Initial phase difference

% Phasing tori 
eps = 1e-5; 
[Str, state(1), Tref] = ICP_guidance(mu, target_orbit.Period, dtheta, k, [r_t0 r_c0], eps, tol);

% ILG transfer 
constraint.Flag = false;
[~, dV, state(2)] = CMLG_guidance(mu, Ln, gamma, target_orbit.Period, constraint, [r_t0 r_c0], tol);

% Periodicity check 
extra = floor(mod(k*Str.Period,target_orbit.Period)/target_orbit.Period*size(target_orbit.Trajectory,1));
chaser_orbit.Trajectory = repmat(chaser_orbit.Trajectory(1:end-1,:), floor(k*Str.Period/target_orbit.Period), 1);
target_orbit.Trajectory = repmat(target_orbit.Trajectory(1:end-1,:), floor(k*Str.Period/target_orbit.Period), 1); 
target_orbit.Trajectory = [target_orbit.Trajectory; target_orbit.Trajectory(1:extra,:)];

% Final rendezvous difference 
S0 = chaser_orbit.Trajectory(1,1:3);                        % Initial guidance point 
Sf = Str.Trajectory(end,1:3)+Str.Trajectory(end,7:9);       % Final guidance point 
e(1) = norm(target_orbit.Trajectory(1,1:3)-S0);             % Initial position space error
e(2) = norm(target_orbit.Trajectory(end,1:3)-Sf);           % Final position space error

% Cost of the maneuver 
dV0 = norm(Str.Trajectory(1,4:6)+Str.Trajectory(1,10:12)-r_c0(4:6));
dVf = norm(target_orbit.Trajectory(end,4:6)-Str.Trajectory(end,4:6)-Str.Trajectory(end,10:12));

%% Results 
%Plot results 
figure(1) 
view(3) 
hold on
T = plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), 'b'); 

% Phasing orbit 
phasing_orbit.Period = Str.Period;
for i = size(Str.Trajectory,1)-1:-1:1
    % Integration of the quasi-periodic trajectory
    tspan = 0:dt:Str.Period;
    [~, St] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, Str.Trajectory(i,:), options);
    S = St(:,1:6)+St(:,7:12);
    guidance = plot3(S(:,1), S(:,2), S(:,3), 'g');
    guidance.Color(4) = 0.2; 
    legend('Target orbit', 'Chaser quasi-periodic orbit', 'AutoUpdate', 'off')
end  

% Initial and final conditions
scatter3(chaser_orbit.Trajectory(1,1), chaser_orbit.Trajectory(1,2),chaser_orbit.Trajectory(1,3), 'filled', 'k')
scatter3(target_orbit.Trajectory(1,1), target_orbit.Trajectory(1,2), target_orbit.Trajectory(1,3), 'filled', 'k'); 
scatter3(Str.Trajectory(end,1)+Str.Trajectory(end,7), Str.Trajectory(end,2)+Str.Trajectory(end,8), Str.Trajectory(end,3)+Str.Trajectory(end,9), 'filled', 'r');
scatter3(target_orbit.Trajectory(end,1), target_orbit.Trajectory(end,2), target_orbit.Trajectory(end,3), 'filled', 'r');
hold off

xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;