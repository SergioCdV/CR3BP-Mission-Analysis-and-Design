%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 17/04/21 % 

%% GNC 5: SMC control law %% 
% This script provides an interface to test the robustness of the continuous control schemes. 

% The relative motion of two spacecraft in the two different halo orbit around L1 and L2 in the
% Earth-Moon system is analyzed, where the chaser is obliged to follow an
% heteroclinic chain between the halos.

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
%Halo characteristics 
Az = 120e6;                                                         %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Homoclinic connection
%Manifold definition
Branch = ['L' 'R'];             %Directions to propagate the manifolds
TOF = tf;                       %Time of flight
rho = 30;                       %Manifold fibers to compute 
 
%Computation flags
long_rendezvous = true;         %Flag to allow for long rendezvous
position_fixed = false;         %Flag to determine a final target state
graphics = false;               %Flag to plot the manifolds

%Final target state
target_orbit.TargetState = shiftdim(target_orbit.Trajectory(1,1:n)); 

%Trajectory design core
[Sg, dV] = HCNC_guidance(mu, Branch, rho, Ln, target_orbit, TOF, long_rendezvous, position_fixed, graphics);;

%% Generate a 2B orbit solution before performing the transfer between halo orbits 
%Orbit conditions
index_1 = Sg.Index;             %Initial position vector epoch
index_2 = Sg.Index+1;           %Final position vector epoc
tof = 1/28;                     %Time of flight for the two body orbit 

s(:,1) = Sg.Trajectory(index_1,1:3).'-[1-mu; 0; 0];      %Initial relative state vector to the Moon 
s(:,2) = Sg.Trajectory(index_2,1:3).'-[1-mu; 0; 0];      %Final relative state vector to the Moon 

%Compute the trajectory arc (this is not being computed accurately)
%[v1, v2, path] = lambert_solver(Mu, s(1:3,1), s(1:3,2));

%Guidance law
Sg.Trajectory = Sg.Trajectory; 

%Integration of the target trajectory
tspan = 0:dt:dt*(size(Sg.Trajectory,1)-1);
r_t0 = target_orbit.Trajectory(1,1:6);                                               %Initial target conditions
[~, St0] = ode113(@(t,s)cr3bp_equations(mu, 1, false, t, s), tspan, r_t0, options);  %Natural target trajectory

%Regression of the guidance coefficients 
Sgr = Sg.Trajectory-St0;                    %Relative guidance law
order = 300; 
[Cp, Cv, Cg, Ci] = CTR_guidance(order, tspan, Sgr);

%Reconstructed guidance trajectory
u = (2*tspan-(tspan(end)+tspan(1)))/(tspan(end)-tspan(1));          %Normalized time domain
T = CH_basis('first', order, u);                                    %Polynomial basis

%Error in the regression
p = Cp*T;                   %Position regression
v = Cv*T;                   %Velocity regression
Sr = St0+[p.' v.'];         %Regress the phase space trajectory

%% GNC algorithms definition 
%Navigation architecture
GNC.Algorithms.Navigation = '';                 %Navigation algorithm

%Control parameters
GNC.Algorithms.Control = 'SMC';                 %Control algorithm

switch (GNC.Algorithms.Control)
    case 'SMC'
        GNC.Algorithms.Solver = 'Encke';                %Dynamics vector field to be solved
        GNC.Guidance.Dimension = 9;                     %Dimension of the guidance law
        GNC.Control.Dimension = 3;                      %Dimension of the control law
        GNC.System.mu = mu;                             %System reduced gravitational parameter
        GNC.Control.SMC.Parameters = [1 1 0.9 1e-3];    %Controller parameters
    case 'LQR'
        GNC.System.mu = mu;                             %System reduced gravitational parameter
        GNC.System.Libration = [Ln gamma];              %Libration point ID
        GNC.Control.LQR.Model = 'RLM';                  %LQR model
        GNC.Control.LQR.Q = 2*eye(9);                   %Penalty on the state error
        GNC.Control.LQR.M = eye(3);                     %Penalty on the control effort
        GNC.Control.LQR.Reference = St0(end,1:3);       %Penalty on the control effort
    case 'SDRE'
        GNC.System.mu = mu;                             %System reduced gravitational parameter
        GNC.System.Libration = [Ln gamma];              %Libration point ID
        GNC.Control.SDRE.Model = 'RLM';                 %SDRE model
        GNC.Control.SDRE.Q = 2*eye(9);                  %Penalty on the state error
        GNC.Control.SDRE.M = eye(3);                    %Penalty on the control effort
    otherwise 
        error('Impulsive guidance tracking has not been implemented yet'); 
end

%Guidance parameters 
GNC.Algorithms.Guidance = 'CTR';               	    %Guidance algorithm
GNC.Guidance.CTR.Order = order;                     %Order of the approximation
GNC.Guidance.CTR.TOF = tspan(end);                  %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	%Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;         %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;     %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.IntergralCoefficients = Ci;        %Coefficients of the Chebyshev approximation

%% GNC: SMC control law for the transfer phase %%
%Initial conditions 
r_c0 = Sg.Trajectory(1,1:6);        %Initial chaser conditions 
rho0 = r_c0-r_t0;                   %Initial relative conditions
s0 = [r_t0 rho0].';                 %Initial conditions of the target and the relative state

switch (GNC.Algorithms.Control)
    case 'LQR'
        s0 = [s0; zeros(3,1)]; 
    case 'SDRE'
        s0 = [s0; zeros(3,1)]; 
end

%Re-integrate trajectory
tic
[t, St] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Control law
[Sgui, ~, u] = GNC_handler(GNC, St(:,1:6), St(:,7:12), t);   

%Error in time 
Stg = [St(:,1:6) St(:,7:end)-Sgui(:,1:6)];
[e, merit] = figures_merit(t, Stg);

%Control effort 
effort = control_effort(t, u, false);

%% GNC: SMC control law for the rendezvous phase %%
%Re-integrate trajectory
GNC.Algorithms.Guidance = '';
tspan2 = 0:dt:2*pi;
tic
[t, St2] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan2, St(end,:), options);
toc

%Complete trajectory
St = [St; St2];
tspan = [tspan tspan(end)+t.'];

%Control law without any guidance law
[~, ~, u2] = GNC_handler(GNC, St2(:,1:6), St2(:,7:12), t);   
u = [u u2];

%Error in time 
[e2, merit2] = figures_merit(t, St2);
e = [e; e2];

%Control effort 
effort2 = control_effort(t, u2, false); 

%% Results %% 
%Plot the transfer
close all; 
figure(1) 
view(3)
hold on 
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), 'b');
plot3(Sg.Trajectory(:,1), Sg.Trajectory(:,2), Sg.Trajectory(:,3), 'k');
plot3(Sr(:,1), Sr(:,2), Sr(:,3), 'r');
hold off
grid on;
title('Heteroclinic guidance orbit between halo orbits')
legend('Target halo orbit', 'Homoclinic orbit', 'Guidance trajectory', 'Location', 'northeast')
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

figure(5)
view(3) 
hold on
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
legend('Rendezvous arc', 'Target orbit', 'Guidance trajectory', 'Location', 'northeast');
title('Converged rendezvous trajectory in the absolute configuration space');

%Configuration space error 
figure(4)
plot(tspan, log(e));
xlabel('Nondimensional epoch');
ylabel('Absolute error $\log{e}$');
grid on;
title('Absolute rendezvous error in the configuration space');

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