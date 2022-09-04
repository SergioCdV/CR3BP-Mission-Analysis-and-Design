%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 24/06/21 % 

%% Potential function study %% 
% This script provides an interface to visualize the relative potential function 
% for periodic target motion. 

% The relative motion of two spacecraft in two halo orbits around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion 

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
tmax = pi;                          %Maximum time of integration (corresponding to a synodic period)
tspan = 0:dt:tmax;                  %Integration time span

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
Az = 200e6;                                                 %Orbit amplitude out of the synodic plane
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Modelling in the synodic frame %% 
%Relative synodic meshgrid
x = -0.1:1e-3:0.1;          %Relative synodic x coordinate
y = x;                      %Relative synodic y coordinate
z = x;                      %Relative synodic z coordinate
[X,Y] = meshgrid(x,y);      %3D meshgrid

%Target orbit motion
rt = target_orbit.Trajectory(:,1:3);

%Evaluate the potential function
U = zeros(size(rt,1), length(x), length(y));     %Preallocation of the potential function

for t = 1:200
    for i = 1:length(x)
        for j = 1:length(y)
                r = [x(i); y(i); 0];                         %Relative position vector
                s = [rt(t,:).'; zeros(3,1); r; zeros(3,1)];  %Augmented phase space vector
                U(t,i,j) = rel_potential(mu, s);             %Potential function
        end
    end
end

%% Results %% 
% Plot results 
rho = 20;
figure(1)
title('Isolines of the potential function'); 
xlabel('Relative synodic x coordinate'); 
ylabel('Relative synodic y coordinate');
for i = 1:size(rt,1)
    contour(X, Y, shiftdim(U(t,:,:)), rho);
    pause(0.1)
end

