%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 14/02/21
% File: test_4.m 
% Issue: 0 
% Validated: 

%% Test 4 %%
% This scripts provides a test interface for the rest of the library
% functions. 

% Test 4 is concerned with validating the single-parameter continuation method, in this case developping 
% a halo family.

%% Test values and constants
%Set graphical environment 
set_graphics(); 

%Numerical setup 
format long; 
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Integration tolerances

%Initial conditions
mu = 0.0121505856;                                          %Reduced gravitational parameter of the system (Earth-Moon)
L = libration_points(mu);                                   %System libration points
Az = 200e6;                                                 %Orbit amplitude out of the synodic plane. Play with it! 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param = [-1 Az Ln gamma m];                                 %Halo orbit parameters (-1 being for southern halo)

%Correction parameters 
dt = 1e-3;                  %Time step to integrate converged trajectories
maxIter = 50;               %Maximum allowed iterations in the differential correction schemes
tol = 1e-10;                %Tolerance 
num = 15;                   %Number of orbits to continuate
direction = -1;             %Direction to continuate (to the Earth)
   
%% Functions
%Compute seed
[halo_seed, haloT] = object_seed(mu, param, 'Halo');        %Generate a halo orbit seed

%Continuation procedure 
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, haloT};                       %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
setup = [mu maxIter tol direction];                         %General setup

[x, state] = continuation(num, method, algorithm, object, corrector, setup);

%% Plotting and results 
%Detect bifurcations 
for i = 1:num
    if (x.Stability(i))
        fprintf('Orbit %i bifurcated. Stability index of %i.\n', i, x.Stability(i));
    end
end

%Plot results
figure(1) 
hold on
tspan = 0:dt:haloT;
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, halo_seed(1,1:6), options);
plot3(S(:,1), S(:,2), S(:,3), 'k');
for i = 1:num
  tspan = 0:dt:x.Period(i);
  [~, S] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, x.Seeds(i,:), options);
  plot3(S(:,1), S(:,2), S(:,3));
end
hold off
xlabel('Synodic normalized x coordinate');
ylabel('Synodic normalized y coordinate');
zlabel('Synodic normalized z coordinate');
title('Converged family of orbits');
grid on;

%% Auxiliary functions 
function set_graphics()
    %Set graphical properties
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
    set(groot, 'defaultAxesFontSize', 11); 
    set(groot, 'defaultAxesGridAlpha', 0.3); 
    set(groot, 'defaultAxesLineWidth', 0.75);
    set(groot, 'defaultAxesXMinorTick', 'on');
    set(groot, 'defaultAxesYMinorTick', 'on');
    set(groot, 'defaultFigureRenderer', 'painters');
    set(groot, 'defaultLegendBox', 'off');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultLegendLocation', 'best');
    set(groot, 'defaultLineLineWidth', 1); 
    set(groot, 'defaultLineMarkerSize', 3);
    set(groot, 'defaultTextInterpreter','latex');
end