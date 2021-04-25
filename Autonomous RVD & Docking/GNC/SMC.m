%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 17/04/21 % 

%% GNC 5: SMC control law %% 
% This script provides an interface to test the Sliding Mode Control strategy for rendezvous missions. 

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% The SMC controller is solved to drive the relative phase space vector to the origin (rendezvous condition).

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

%Spacecraft mass 
mass = 1e-10;

%Time span 
dt = 1e-3;                          %Time step
tf = 2*pi;                          %Rendezvous time
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
Az = 200e6;                                                 %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                    %Number of orbits to continuate
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};         %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
direction = 1;                                              %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                         %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(2,:), maxIter, tol);

%% Modelling in the synodic frame %%
index = fix(tf/dt);                                         %Rendezvous point
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% GNC: discrete SMC control law %%
%Preallocation 
Sc = zeros(length(tspan), 2*n);                             %Relative orbit trajectory
Sc(1,:) = Sn(1,:);                                          %Initial relative state
u = zeros(3,length(tspan));                                 %Control law
e = zeros(1,length(tspan));                                 %Error to rendezvous 

%Reference state 
refState = zeros(n+3,1);                                    %Reference state (rendezvous condition)

%Compute the trajectory
for i = 1:length(tspan)
    %Compute the control law 
    u(:,i) = hybrid_controller(mu, refState, Sc(i,:));

    %Re-integrate trajectory
    [~, s] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke C', t, s, u(:,i)), [0 dt], S(i,:), options);
    
    %Update initial conditions
    Sc(i+1,:) = s(end,:);
        
    %Error in time 
    e(i) = norm(s(end,7:12));
end

%% GNC: SMC control law %%
% %Preallocation 
e = zeros(1,length(tspan));                                 %Error to rendezvous 

%Reference state 
refState = zeros(n+3,1);                                    %Reference state (rendezvous condition)

%Re-integrate trajectory
[~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke SMC', t, s, refState), tspan, s0, options);

%Error in time 
for i = 1:length(tspan)
    e(i) = norm(Sc(i,7:12));
end
    
%% Results %% 
%Plot results 
figure(1) 
view(3) 
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3)); 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
hold off
legend('Target motion', 'Chaser motion'); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Reconstruction of the natural chaser motion');

%Plot relative phase trajectory
figure(2) 
view(3) 
plot3(Sc(:,7), Sc(:,8), Sc(:,9)); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative motion in the configuration space');

%Configuration space error 
figure(3)
plot(tspan, e); 
xlabel('Nondimensional epoch');
ylabel('Absolute error');
grid on;
title('Absolute error in the configuration space (L2 norm)');

%Rendezvous animation 
if (true)
    figure(4) 
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

%% Auxiliary functions 
%Control algorithms
function [T] = hybrid_controller(mu, refState, x)
    %SMC parameters 
    lambda = 1;                %General loop gain
    epsi = 1;                  %Reachability condition gain
    alpha = 0.9;               %Reachability condition exponent
    delta = 1e-2;              %Boundary layer width
    
    %Attitude state
    r = x(7:9).';              %Instanteneous position vector
    v = x(10:12).';            %Instanteneous velocity vector
    
    %Compute the position and velocity errors
    dr = r-refState(1:3);      %Position error
    dv = v-refState(4:6);      %Velocity error
    
    %Torque computation
    s = dv+lambda*dr;                                                   %Sliding surface
    f = nlr_model(mu, true, false, 'Encke C', 0, x.', zeros(3,1));      %Relative CR3BP dynamics
    ds = epsi*(norm(s)^(alpha)*saturation(s, delta).'+s);               %Reachability condition function
    T = refState(7:9)-f(10:12)-lambda*dv-ds;                            %Control vector
end

%Saturation function
function [U] = saturation(s, delta)
    %Compute a bang bang saturation law, given the boundary layer delta 
    U = zeros(1,length(s));
    for i = 1:size(s)
        if (s(i) > delta)
            U(i) = 1; 
        elseif (s(i) < -delta)
            U(i) = -1; 
        else
            U(i) = (1/delta)*s(i);
        end
    end
end