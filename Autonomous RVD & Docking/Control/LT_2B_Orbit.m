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

%Optimization 
optimization = false;                %Optimize the controller parameters 

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
%Target halo characteristics 
Az = 120e6;                                                         %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, ~] = object_seed(mu, halo_param, 'Halo');               %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Target halo characteristics 
Az = 120e6;                                                         %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 2;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, ~] = object_seed(mu, halo_param, 'Halo');               %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[initial_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Heteroclinic connection
%Manifold definition
TOF = 2*pi;                     %Time of flight
rho = 30;                       %Manifold fibers to compute 

%Computation flags
position_fixed = false;         %Flag to determine a final target state
graphics = true;                %Flag to plot the manifolds

%Final target state
target_orbit.TargetState = shiftdim(target_orbit.Trajectory(1,1:n)); 

%Connection itenerary
sequence = [1 2];       %Connection itenerary
branch = ['R' 'L'];     %Manifold branches to be propagated

%Trajectory design core
[Sg, dV] = HTRC_guidance(mu, sequence, branch, rho, target_orbit, initial_orbit, TOF, position_fixed, graphics);

%% Generate a 2B orbit solution before performing the transfer between halo orbits 
%Orbit conditions
r_m2 = Sg.Trajectory(Sg.Index,1:3).'-[1-mu; 0; 0];      %Relative position vector to the Moon 
r_m2 = norm(r_m2);                                      %Relative orbit semimajor axis
theta = 0:1e-2:2*pi;                                    %Two body anomaly 
v_m2 = sqrt(mu/r_m2);                                   %Velocity of the circular orbit 

%Guidance law
Sg.Trajectory = [Sg.Trajectory(1:Sg.Index,:); zeros(length(theta),6); Sg.Trajectory(Sg.Index+1:end,:)]; 

for i = 1:length(theta)
    %Position vector
    Sg.Trajectory(Sg.Index+i,1:3) = [1-mu; 0; 0].' + norm(r_m2)*[cos(theta(i)) sin(theta(i)) 0]; 

    %Velocity vector 
    Sg.Trajectory(Sg.Index+i,4:6) = v_m2*[-sin(theta(i)) cos(theta(i)) 0];                        
end

%Integration of the target trajectory
tspan = 0:dt:dt*(size(Sg.Trajectory,1)-1);
r_t0 = target_orbit.Trajectory(1,1:6);                                               %Initial target conditions
[~, St0] = ode113(@(t,s)cr3bp_equations(mu, 1, false, t, s), tspan, r_t0, options);  %Natural target trajectory

%Regression of the guidance coefficients 
Sgr = Sg.Trajectory-St0;                    %Relative guidance law
order = 300; 
[Cp, Cv, Cg] = CTR_guidance(order, tspan, Sgr);

%Reconstructed guidance trajectory
T = zeros(order, length(tspan));                                    %Preallocation of the polynomial basis
u = (2*tspan-(tspan(end)+tspan(1)))/(tspan(end)-tspan(1));          %Normalized time domain

for i = 1:length(tspan)
    T(:,i) = chebyshev('first', order, u(i));
end

%Error in the regression
p = Cp*T;                   %Position regression
v = Cv*T;                   %Velocity regression
Sr = St0+[p.' v.'];         %Regress the phase space trajectory

%% GNC algorithms definition 
%Control parameters
GNC.Algorithms.Navigation = '';                 %Navigation algorithm
GNC.Algorithms.Control = 'SMC';                 %Control algorithm
GNC.Algorithms.Solver = 'Encke';                %Dynamics vector field to be solved
GNC.Guidance.Dimension = 9;                     %Dimension of the guidance law
GNC.Control.Dimension = 3;                      %Dimension of the control law
GNC.System.mu = mu;                             %System reduced gravitational parameter
GNC.Control.SMC.Parameters = [100 1 0.9 1e-3];     %Controller parameters

%Guidance parameters 
GNC.Algorithms.Guidance = 'CTR';               	    %Guidance algorithm
GNC.Guidance.CTR.Order = order;                     %Order of the approximation
GNC.Guidance.CTR.TOF = tspan(end);                  %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	%Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;         %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;     %Coefficients of the Chebyshev approximati

%% GNC: SMC control law for the transfer phase %%
%Initial conditions 
r_c0 = initial_orbit.Trajectory(1,1:6);            %Initial chaser conditions 
rho0 = r_c0-r_t0;                                  %Initial relative conditions
s0 = [r_t0 rho0].';                                %Initial conditions of the target and the relative state

%Re-integrate trajectory
tic
[t, St] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Control law
[Sgui, ~, u] = GNC_handler(GNC, St(:,1:6), St(:,7:12), t);   

%Error in time 
Stg = [St(:,1:6) St(:,7:end)-Sgui(:,1:6)];
[e, merit] = figures_merit(tspan, Stg);

%Control effort 
effort = control_effort(t, u, false);

%% GNC: SMC control law for the rendezvous phase %%
%Re-integrate trajectory
GNC.Algorithms.Guidance = '';
tic
[t, St2] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, St(end,:), options);
toc

%Complete trajectory
St = [St; St2];
tspan = [tspan tspan(end)+tspan];

%Control law without any guidance law
[~, ~, u] = GNC_handler(GNC, St2(:,1:6), St2(:,7:12), t);   

%Error in time 
[e2, merit2] = figures_merit(t, St2);
e = [e; e2];

%Control effort 
effort2 = control_effort(t, u, false); 

%% Results %% 
%Plot the transfer
close all; 
figure(1) 
view(3)
hold on 
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), 'b');
plot3(initial_orbit.Trajectory(:,1), initial_orbit.Trajectory(:,2), initial_orbit.Trajectory(:,3), 'b');
plot3(Sg.Trajectory(:,1), Sg.Trajectory(:,2), Sg.Trajectory(:,3), 'k');
plot3(Sr(:,1), Sr(:,2), Sr(:,3), 'r');
hold off
grid on;
title('Heteroclinic guidance orbit between halo orbits')
legend('Target halo orbit', 'Initial halo orbit', 'Heteroclinic orbit', 'Guidance trajectory', 'Location', 'northeast')
xlabel('Synodic $x$ coordinate')
ylabel('Synodic $y$ coordinate')
zlabel('Synodic $z$ coordinate')

%Plot relative phase trajectory
figure(2) 
view(3) 
plot3(St(:,7), St(:,8), St(:,9)); 
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
title('Motion in the relative configuration space');

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan, St(:,7)); 
plot(tspan, St(:,8)); 
plot(tspan, St(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinates');
grid on;
legend('$x$', '$y$', '$z$');
title('Relative position in time');
subplot(1,2,2)
hold on
plot(tspan, St(:,10)); 
plot(tspan, St(:,11)); 
plot(tspan, St(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinates');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Relative velocity in time');

%Configuration space error 
figure(4)
plot(tspan, log(e));
xlabel('Nondimensional epoch');
ylabel('Absolute error $\log{e}$');
grid on;
title('Absolute rendezvous error in the configuration space');

figure(5)
view(3) 
hold on
c = plot3(initial_orbit.Trajectory(:,1), initial_orbit.Trajectory(:,2), initial_orbit.Trajectory(:,3), 'r', 'Linewidth', 0.1); 
r = plot3(St(:,7)+St(:,1), St(:,8)+St(:,2), St(:,9)+St(:,3), 'b', 'Linewidth', 0.1); 
t = plot3(St0(:,1), St0(:,2), St0(:,3), 'r', 'Linewidth', 0.1);
g = plot3(Sg.Trajectory(:,1), Sg.Trajectory(:,2), Sg.Trajectory(:,3), 'k', 'Linewidth', 0.1);
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
hold off
text(L(1,Ln)+1e-3, L(2,Ln), 0, '$L_1$');
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Initial orbit', 'Rendezvous arc', 'Target orbit', 'Guidance trajectory', 'Location', 'northeast');
title('Converged rendezvous trajectory in the absolute configuration space');

%Rendezvous animation 
if (false)
    figure(5) 
    view(3) 
    grid on;
    hold on
    plot3(initial_orbit.Trajectory(:,1), initial_orbit.Trajectory(:,2), initial_orbit.Trajectory(:,3), 'k-.'); 
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