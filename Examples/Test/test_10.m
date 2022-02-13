%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/10/21
% File: test_10.m 
% Issue: 0 
% Validated: 

%% Test 10 %%
% This scripts provides a test interface for the rest of the library
% functions. 

% Test 10 is concerned with the Encke's method integration of periodic orbits.

%% Test values and constants
%Set graphical environment 
set_graphics(); 

%Initial conditions
mu = 0.0121505856;                                          %Reduced gravitational parameter of the system (Earth-Moon)
L = libration_points(mu);                                   %System libration points
Az = 60e1;                                                 %Orbit amplitude out of the synodic plane. Play with it!
Ax = 60e1;                                                 %Orbit amplitude in the synodic plane. Play with it! 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ax = dimensionalizer(384400e3, 1, 1, Ax, 'Position', 0);    %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param_halo = [1 Az Ln gamma m];                             %Halo orbit parameters (-1 for southern halo)

%Correction parameters 
maxIter = 50;      %Maximum allowed iterations in the differential correction schemes
tol = 1e-10;       %Tolerance 

%Numerical setup 
format long; 
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Integration tolerances

%% Functions
%Compute seeds
[halo_seed, haloT] = object_seed(mu, param_halo, 'Halo');   %Generate a halo orbit seed

%Halo orbit 
[halo_orbit, state(2)] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Initial conditions and time span 
K = 5;                                   %Number of periods to integrate
s0 = halo_orbit.Trajectory(1,1:6).';     %Initial conditions
dt = 1e-3;                               %Time step
tspan = 0:dt:K*halo_orbit.Period;        %Integration time span

%Propagation setup 
setup.Method = 'Newton';

%Newton integration
tic
[~, S_N] = ode113(@(t,s)cr3bp_propagator(setup, mu, true, false, t, s), tspan, s0, options);
toc

%Encke's integration 
setup.Method = 'Encke';
S0 = s0; 

S_E = zeros(size(S_N));

index = [1 2 4 5];
for i = 1:length(index)
    setup.LagrangePointID = index(i);
    setup.LagrangePointPosition = L(1:3,setup.LagrangePointID);
    
    s0 = S0; 
    s0(1:3) = s0(1:3)-L(1:3,setup.LagrangePointID);
    tic
    [t, Saux] = ode113(@(t,s)cr3bp_propagator(setup, mu, true, false, t, s), tspan, s0, options);
    toc
    S_E = S_E + Saux; 
    S_E(:,1:3) = S_E(:,1:3)+L(1:3,setup.LagrangePointID).'; 
end
S_E = S_E/i;

%Compute the error with respect to the differentially corrected orbit 
e = zeros(2,size(S_E,1));       %Preallocation 

k = 1; 
for i = 1:size(S_E,1)
    if (mod(t(i),halo_orbit.Period) == 0)
        k = 1;
    else
        k = k+1;
    end

    e(1,i) = norm(halo_orbit.Trajectory(k,1:3)-S_N(i,1:3))/norm(halo_orbit.Trajectory(k,1:3));
    e(2,i) = norm(halo_orbit.Trajectory(k,1:3)-S_E(i,1:3))/norm(halo_orbit.Trajectory(k,1:3));
end

%% Plotting
%Plot the three orbits
figure(1) 
view(3)
hold on
plot3(halo_orbit.Trajectory(:,1), halo_orbit.Trajectory(:,2), halo_orbit.Trajectory(:,3), 'g');
plot3(S_N(:,1), S_N(:,2), S_N(:,3), 'b')
plot3(S_E(:,1), S_E(:,2), S_E(:,3), 'r')
hold off
xlabel('Synodic normalized $x$ coordinate');
ylabel('Synodic normalized $y$ coordinate');
zlabel('Synodic normalized $z$ coordinate');
title('Numerically integrated halo orbits');
legend('Differential corrected', 'Newton', 'Encke');
grid on;

%Plot the error 
figure(2) 
hold on 
plot(tspan, log(e(1,:)), 'b'); 
plot(tspan, log(e(2,:)), 'r')
hold off
grid on; 
title('Relative error to the exact solution');
legend('Newtonian error', 'Encke error');

