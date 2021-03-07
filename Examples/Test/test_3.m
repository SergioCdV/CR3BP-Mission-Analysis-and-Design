%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/02/21
% File: test_3.m 
% Issue: 0 
% Validated: 

%% Test 3 %%
% This scripts provides a test interface for the rest of the library
% functions. 

% Test 3 is concerned with the computation of invariant manifols associated
% with periodic orbits.

%% General setup 
format long;                                                %Set command window numerical format
set_graphics();                                             %Set up the graphical environment

%% Initial conditions 
mu = 0.0121505856;                                          %Reduced gravitational parameter of the system
L = libration_points(mu);                                   %Libration points computation

%Generate a Lyapunov orbit seed
Az = 50e6;                                                  %Orbit amplitude out of the synodic plane. Play with it!
Ax = 50e6;                                                  %Orbit amplitude in the synodic plane. Play with it! 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ax = dimensionalizer(384400e3, 1, 1, Ax, 'Position', 0);    %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param = [Ax Az 0 0 Ln gamma m];                             %Lyapunov orbit parameters
lyapunov_seed = object_seed(mu, param, 'Lyapunov');         %Lyapunov seed     

%Time integration
Tf = 20*pi;                                                 %Final time
dt = 0.01;                                                  %Time step of 1 second
tspan = 0:dt:Tf;                                            %Integration time span

%Differential correction scheme set up
maxIter = 50;                                               %Maximum allowed iterations in the differential correction schemes
tol = 1e-5;                                                 %Tolerance 

%% Main test
%Orbit computation
[lyapunov_orbit, state] = differential_correction('Planar', mu, lyapunov_seed, maxIter, tol);

%Manifold computation
rho = 40;                    %Number of manifold fibers to compute
manifold_ID = 'U';           %Unstable manifold (U or S)
manifold_branch = 'L';       %Left branch of the manifold (L or R)

Manifold = invariant_manifold(mu, manifold_ID, manifold_branch, lyapunov_orbit.Trajectory, rho, 0:dt:2*pi);

%% Plotting and results 
figure(1) 
view(3)
plot3(lyapunov_orbit.Trajectory(:,1), lyapunov_orbit.Trajectory(:,2), lyapunov_orbit.Trajectory(:,3));
xlabel('Synodic normalized x coordinate');
ylabel('Synodic normalized y coordinate');
zlabel('Synodic normalized z coordinate');
title('Converged Lyapunov orbit');
grid on;

figure(2)
hold on 
for i = 1:size(Manifold,1)
    ManifoldAux = shiftdim(Manifold(i,:,:));
    plot3(ManifoldAux(:,1), ManifoldAux(:,2), ManifoldAux(:,3), 'm');
end
plot3(lyapunov_orbit.Trajectory(:,1), lyapunov_orbit.Trajectory(:,2), lyapunov_orbit.Trajectory(:,3));
% plot3(L(1,1), L(2,1), L(3,1), 'ok');
% plot3(L(1,2), L(2,2), L(3,2), 'ok');
% plot3(L(1,3), L(2,3), L(3,3), 'ok');
% plot3(L(1,4), L(2,4), L(3,4), 'ok');
% plot3(L(1,5), L(2,5), L(3,5), 'ok');
% plot3(-mu, 0, 0, 'ob');
% plot3(1-mu, 0, 0, 'ob');
hold off
xlabel('Synodic normalized x coordinate');
ylabel('Synodic normalized y coordinate');
zlabel('Synodic normalized z coordinate');
title('Lyapunov orbit and associated manifold');
grid on;
