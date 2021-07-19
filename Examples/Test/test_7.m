%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/07/21
% File: test_7.m 
% Issue: 0 
% Validated: 

%% Test 7 %%
% This scripts provides a test interface for the rest of the library
% functions. 

% Test 7 is concerned with differential correction algorithms for computation of dynamical 2D tori.

%% Test values and constants
%Set graphical environment 
set_graphics(); 

%Initial conditions
mu = 0.0121505856;                                          %Reduced gravitational parameter of the system (Earth-Moon)
L = libration_points(mu);                                   %System libration points
Az = 50e6;                                                  %Orbit amplitude out of the synodic plane. Play with it!
Ax = 50e6;                                                  %Orbit amplitude in the synodic plane. Play with it! 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ax = dimensionalizer(384400e3, 1, 1, Ax, 'Position', 0);    %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param_halo = [1 Az Ln gamma m];                             %Halo orbit parameters (-1 for southern halo)
param_lyap = [Ax Az 0 0 Ln gamma m];                        %Lyapunov orbit parameters

%Correction parameters 
maxIter = 20;      %Maximum allowed iterations in the differential correction schemes
tol = 1e-4;        %Tolerance 

%% Functions
%Compute seeds
[halo_seed, haloT] = object_seed(mu, param_halo, 'Halo');   %Generate a halo orbit seed
lyapunov_seed = object_seed(mu, param_lyap, 'Lyapunov');    %Generate a Lyapunov orbit seed
vertical_seed = [0.9261 0 0.3616 0 -0.0544  0];             %State vector of a vertical orbit

%Generate the orbits
[lyapunov_orbit, ~] = differential_correction('Planar', mu, lyapunov_seed, maxIter, tol);
[halo_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
[vertical_orbit, ~] = differential_correction('Double Plane Symmetric', mu, vertical_seed, maxIter, tol);

%Compute the halo tori 
J = jacobi_constant(mu, halo_orbit.Trajectory(end,1:6).');
[xf, state] = differential_torus('Single shooting energy', mu, halo_orbit, maxIter, tol, 17, J);

%% Surface interpolation 
[theta, phi] = meshgrid(0:1e-2:2*pi,0:1e-2:2*pi);
x = shiftdim(xf.Trajectory(:,:,1));
y = shiftdim(xf.Trajectory(:,:,2));
z = shiftdim(xf.Trajectory(:,:,3));

%% Plot results 
figure(1) 
view(3)
hold on 
for i = 1:size(xf.Trajectory,1)
    S = shiftdim(xf.Trajectory(i,:,:));
    plot3(S(:,1), S(:,2), S(:,3));
end
hold off
grid on;
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
title('Computed torus');