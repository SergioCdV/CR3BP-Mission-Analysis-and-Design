%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 22/12/21 % 

%% ISTS Gateway resupply mission %% 
% This script provides an interface to test the robustness of the continuous control schemes. 

% A resupply dynamical chain is established for nearly continuous supply to an L2 spacecraft.

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
optimization = false;               %Optimize the controller parameters 

%Time span 
dt = 1e-3;                          %Time step
tf = 2*pi;                          %Rendezvous time
tspan = 0:dt:tf;                    %Integration time span

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
Ln = 2;                                                             %Orbits around L2
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters);

%Compute the NRHO
[halo_seed, haloT] = object_seed(mu, halo_param, 'Halo');           %Generate a halo orbit seed

%Continuation procedure 
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                    %Number of orbits to continuate
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, haloT};                       %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
direction = 1;                                              %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                         %General setup

[target_orbit, state_energy] = continuation(num, method, algorithm, object, corrector, setup);

%Generate the NRHO
s0 = [target_orbit.Seeds(end,:).'; reshape(eye(n), [n^2 1])];
tspan = 0:dt:target_orbit.Period(end);
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);

%Target Gateway orbit
Sn = S; 

%Loitering halo characteristics 
Az = 120e6;                                                         %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, ~] = object_seed(mu, halo_param, 'Halo');               %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[loiter_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Modelling in the synodic frame %%
%Integration of the model
tspan = 0:dt:target_orbit.Period(end);
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);
Sn = repmat(S, floor(tf/target_orbit.Period(end))+1, 1);     

%Guidance law coefficients
Sg = jacobi_constant(mu,Sn(1,1:n).');

%% GNC algorithms definition 
GNC.Algorithms.Guidance = '';                       %Guidance algorithm
GNC.Algorithms.Navigation = '';                     %Navigation algorithm
GNC.Algorithms.Control = 'MLQR';                    %Control algorithm

GNC.Guidance.Dimension = 9;                         %Dimension of the guidance law
GNC.Control.Dimension = 3;                          %Dimension of the control law

GNC.System.mu = mu;                                 %Systems's reduced gravitational parameter
GNC.Control.MLQR.Q = eye(2);                        %Penalty on the state error
GNC.Control.MLQR.M = eye(3);                        %Penalty on the control effort
GNC.Control.MLQR.Reference = Sg;                    %Penalty on the control effort
GNC.Control.MLQR.Period = target_orbit.Period(end); %Penalty on the control effort

%Floquet mode matrix
[E, lambda] = eig(reshape(Sn(end,n+1:end), [n n]));
GNC.Control.MLQR.FloquetModes = diag(log(diag(lambda))/target_orbit.Period(end));
for i = 1:size(E,2)
    E(:,i) = E(:,i)/lambda(i,i);
end
GNC.Control.MLQR.FloquetDirections = E; 

%% Gateway stationkeeping orbit
%Initial conditions 
k = dimensionalizer(Lem, 1, 1, 10, 'Position', 0);      %Noise gain
r_t0 = target_orbit.Seeds(end,1:n);                     %Initial guidance target conditions
s0 = r_t0-2*k*ones(1,n)+k*rand(1,n);                    %Noisy initial conditions
tspan = 0:dt:tf;                                        %New integration time span

GNC.Navigation.NoiseVariance = 0;                       %Noise variance

%Compute the trajectory
s0 = [r_t0 s0-r_t0 reshape(eye(n), [1 n^2])];
tic
[~, Staux1] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan(tspan < 0.1), s0, options);
[~, Staux2] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan(size(Staux1,1):end), Staux1(end,:), options);
toc
Sgate = [Staux1; Staux2(2:end,:)];

%Control law
[~, ~, u] = GNCc_handler(GNC, Staux1(:,1:n), Staux1(:,n+1:end), tspan(tspan < 0.1));

%Control integrals
energy_sk = control_effort(tspan(tspan < 0.1), u, false);

%% Heteroclinic connection between the loitering and target orbit
%Manifold definition
TOF = 2*pi;                     %Time of flight
rho = 90;                       %Manifold fibers to compute 

%Computation flags
position_fixed = false;         %Flag to determine a final target state
graphics = true;               %Flag to plot the manifolds

%Final target state
target_orbit.Trajectory = Sn;
target_orbit.TargetState = shiftdim(target_orbit.Seeds(end,1:n)); 

%Connection itenerary
sequence = [1 2];               %Connection itenerary
branch = ['L' 'L'];             %Manifold branches to be propagated

%Trajectory design core
[Shr, dVhr] = HTRC_guidance(mu, sequence, branch, rho, target_orbit, loiter_orbit, TOF, position_fixed, graphics);

%Control integrals
energy_hr = control_effort(NaN, dVhr, true);

%Guidance law
Sg = Shr.Trajectory; 

%Integration of the target trajectory
tspan = 0:dt:dt*(size(Sg,1)-1);
r_t0 = target_orbit.Trajectory(1,1:6);                                               %Initial target conditions
[~, St0] = ode113(@(t,s)cr3bp_equations(mu, 1, false, t, s), tspan, r_t0, options);  %Natural target trajectory

%Regression of the guidance coefficients 
Sgr = Sg-St0;                                            %Relative guidance law
order = 300;                                             %Order of the approximation
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

figure(4)
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

%Configuration space error 
figure(5)
plot(tspan, log(e));
xlabel('Nondimensional epoch');
ylabel('Absolute error $\log{e}$');
grid on;
title('Absolute rendezvous error in the configuration space');

%Rendezvous animation 
if (false)
    figure(6) 
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