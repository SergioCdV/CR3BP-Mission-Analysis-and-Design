%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 12: APF guidance-control law %% 
% This script provides an interface to test Artificial Potential Functions guidance laws for
% rendezvous missions.

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
Az = 200e6;                                                         %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                                     %Bifucartion tolerance on the stability index
num = 2;                                                            %Number of orbits to continuate
method = 'SPC';                                                     %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                        %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};                 %Object and characteristics to continuate
corrector = 'Plane Symmetric';                                      %Differential corrector method
direction = 1;                                                      %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                                 %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(2,:), maxIter, tol);

%% Modelling in the synodic frame %%
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                       %Reconstructed chaser motion via Encke method

%% APF guidance scheme
%Obstacles 
So = 10*ones(3,1);                               %Relative positon of the obstalces

%Controller scheme penalties
Penalties.AttractivePenalty = eye(3);            %Penalty on the distance to the origin
Penalties.RepulsivePenalty = eye(3);             %Penalty on the distance to the obstacles 
Penalties.RepulsiveWidth = 1e-3;                 %Width of the repulsive function

%Safety message 
safe_corridor.Safety = false;
safe_corridor.Parameters(1) = deg2rad(10);       %Safety corridor angle
safe_corridor.Parameters(2) = 0;                 %Safety distance to the docking port
safe_corridor.Parameters(3:4) = [1 1];           %Dimensions of the safety corridor

%Compute the guidance law
output = APF_guidance(safe_corridor, Penalties, So, tf, s0(7:12), true);
Cp = output(1:3,:);
Cv = output(4:6,:);
Cg = output(7:9,:);
 
%% GNC trajectory: APF-SMC
%GNC control structure
GNC.Algorithms.Guidance = 'APF';                    %Guidance algorithm
GNC.Algorithms.Navigation = '';                     %Navigation algorithm
GNC.Algorithms.Control = 'SMC';                     %Control algorithm
GNC.Guidance.Dimension = 9;                         %Dimension of the guidance law
GNC.Control.Dimension = 3;                          %Dimension of the control law

GNC.Guidance.CTR.Order = size(Cp,2);                %Order of the approximation
GNC.Guidance.CTR.TOF = tf;                          %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;         %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;         %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;     %Coefficients of the Chebyshev approximation

GNC.Guidance.APF.Safety = safe_corridor;            %APF safety parameters
GNC.Guidance.APF.Penalties = Penalties;             %APF guidance core parameters
GNC.Guidance.APF.Obstacles = So;                    %Relative position of the obstacles
            
GNC.System.mu = mu;                                 %System reduced gravitational parameter
GNC.Control.SMC.Parameters = [10 10 0.9 0.01];      %Controller parameters

%Re-integrate trajectory
[t, St] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);

%Error in time 
[e, merit] = figures_merit(tspan, St);

%Control law
[Sg, ~, u] = GNC_handler(GNC, St(:,1:6), St(:,7:12), t);    

%Control effort 
effort = control_effort(tspan, u, false);

%% Results 
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
hold on
plot3(St(:,7), St(:,8), St(:,9)); 
plot3(Sg(:,1), Sg(:,2), Sg(:,3));
hold off
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative motion in the configuration space');
legend('Controlled trajectory', 'Reference trajectory')

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan, St(:,7)); 
plot(tspan, St(:,8)); 
plot(tspan, St(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('x coordinate', 'y coordinate', 'z coordinate');
title('Relative position evolution');
subplot(1,2,2)
hold on
plot(tspan, St(:,10)); 
plot(tspan, St(:,11)); 
plot(tspan, St(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinate');
grid on;
legend('x velocity', 'y velocity', 'z velocity');
title('Relative velocity evolution');

%Configuration space error 
figure(4)
plot(tspan, log(e)); 
xlabel('Nondimensional epoch');
ylabel('Absolute error  (log)');
grid on;
title('Absolute error in the configuration space (L2 norm)');

%Rendezvous animation 
if (false)
    figure(5) 
    view(3) 
    grid on;
    hold on
    plot3(St(:,1), St(:,2), St(:,3), 'k-.'); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    title('Rendezvous simulation');
    for i = 1:size(St,1)
        T = scatter3(St(i,1), St(i,2), St(i,3), 30, 'b'); 
        V = scatter3(St(i,1)+St(i,7), St(i,2)+St(i,8), St(i,3)+St(i,9), 30, 'r');

        drawnow;
        delete(T); 
        delete(V);
    end
    hold off
end