%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 21/10/22 % 

%% James Webb Mission for MPDI %% 
% This script provides an interface to test the general control scheme for a rendezvous, docking and undocking mission for a JWST mission.

% Units are non-dimensional and solutions are expressed in the synodic
% reference frame as defined by Howell, 1984.

%% Set up %%
clc 
clear; 
close all; 

% Set up graphics 
set_graphics();

% Integration options
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 

%% Contants and initial data %% 
% Phase space dimension 
n = 6; 

% CR3BP constants 
mu = 3.003e-6;                      % Sun-Earth reduced gravitational parameter
Lem = 149597870700;                 % Mean distance from the Sun to the Earth
Tc = 86400*365;                     % Sun-Earth characteristic period
Vc = 29.784e3;                      % Sun-Earth characteristic velocity

% General setup 
dt = 1e-3;                          % Time step to integrate converged trajectories
maxIter = 20;                       % Maximum allowed iterations in the differential correction schemes
tol = 1e-10;                        % Differential correction tolerance 
 
%% Target halo orbit computation %%
% Initial conditions
L = libration_points(mu);                                   % System libration points
Az = 220e6;                                                 % Orbit amplitude out of the synodic plane. Play with it! 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         % Normalize distances for the E-M system
Ln = 2;                                                     % Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          % Li distance to the second primary
m = 1;                                                      % Number of periods to compute
param = [1 Az Ln gamma m];                                  % Halo orbit parameters (-1 being for southern halo)
   
% Compute the nominal orbit
[halo_seed, haloT] = object_seed(mu, param, 'Halo');        % Generate the halo orbit seed

% Generate the orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

% Initial conditions
L = libration_points(mu);                                   % System libration points
Az = 90e6;                                                 % Orbit amplitude out of the synodic plane. Play with it! 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         % Normalize distances for the E-M system
Ln = 2;                                                     % Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          % Li distance to the second primary
m = 1;                                                      % Number of periods to compute
param = [1 Az Ln gamma m];                                  % Halo orbit parameters (-1 being for southern halo)
   
% Compute the nominal orbit
[halo_seed, haloT] = object_seed(mu, param, 'Halo');        % Generate the halo orbit seed

% Generate the orbit
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Natural spacecraft motion %%
r_t0 = target_orbit.Trajectory(500, 1:6);         % Initial target conditions
r_c0 = chaser_orbit.Trajectory(1, 1:6);           % Initial chaser conditions 
rho0 = r_c0-r_t0;                                 % Initial relative conditions
s0 = [r_t0 rho0].';                               % Initial conditions of the target and the relative state

% Integration of the model
tspan = 0:dt:2*pi; 
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
Sn = S;                

% Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                        % Reconstructed chaser motion via Encke method

% Insertion plot 
figure(1)
view(3)
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3));
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3));
scatter3(Sn(1,1), Sn(1,2), Sn(1,3), 'blue', 'filled');
scatter3(S_rc(1,1), S_rc(1,2), S_rc(1,3), 'red', 'filled');
hold off
grid on;
legend('Target orbit', 'Loiter orbit', 'Target at HOI', 'Chaser at HOI'); 
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');

%% Initial phase: MPC + discrete iLQR 
% Algorithm setup
GNC.Control.iLQR.Mode = 'Discrete';                                 % iLQR solver
GNC.Control.LQR.Q = 1e1*[eye(3) zeros(3); zeros(3) 1e-6*eye(3)];        % Penalty on the state error
GNC.Control.LQR.M = eye(3);                                         % Penalty on the control effort

dVmax = 3e-3;                                                       % Maximum control authority

tol = [1e-4 1e-5];                                                  % Convergence tolerance

N = 60;                                 % Number of impulses
dt = 5e-2;                              % Time step for the phase 
tf(1) = (N+1)*dt;                       % Final TOF 
tspan = 0:dt:tf(1);                     % Integration span 

% Preallocation for speed
u = zeros(3,length(tspan)-1);           % Final control sequence
S1 = zeros(length(tspan), 2*n);         % Final rendezvous trajectory

% Noise constants 
Kp = 5/Lem;             % Position uncertainty
Kv = 0.5/Vc;            % Velocity uncertainty

S1(1,:) = s0; 

% Controller
tic
for i = 1:length(tspan)-1
    % Guidance 
    tf_aux = tf(1)-dt*(i-1);
    s0 = S1(i,:);
    s0(n+1:2*n) = s0(n+1:2*n)+[Kp Kp Kp Kv Kv Kv].*normrnd(0,1,1,n);
    [~, St_ilqr, u_aux, ~] = iLQR_control(mu, tf_aux, s0, GNC, dVmax, dt, false, tol);

    % Integration 
    u(:,i) = u_aux(:,1); 

    % Noisy integration
    du = rand(1,3);
    du = 0.05*rand*norm(u(:,1))*du/norm(du);
    U = [zeros(1,n+3) u_aux(:,1).'+du];
    s0 = S1(i,:)+U; 
    [~, aux] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), [0 dt], s0, options);
    S1(i+1,:) = aux(end,:);
end
Time(1) = toc; 

% Error in time 
[e_ilqr, merit_ilqr] = figures_merit(tspan, S1(:,7:12));

% Control integrals
effort_ilqr = control_effort(tspan, u, true);

% Range reduction 
range(1) = 1-norm(S1(end,7:9))/norm(S1(1,7:9));

%% Second phase: SMC-offline LQR guidance 
% Guidance computation
GNC.Algorithms.Guidance = '';                   % Guidance algorithm
GNC.Algorithms.Navigation = '';                 % Navigation algorithm
GNC.Algorithms.Control = 'LQR';                 % Control algorithm

GNC.Guidance.Dimension = 9;                     % Dimension of the guidance law
GNC.Control.Dimension = 3;                      % Dimension of the control law
GNC.Navigation.NoiseVariance = 0;               % Noise variance

GNC.System.mu = mu;                             % Systems's reduced gravitational parameter
GNC.System.Libration = [Ln gamma];              % Libration point ID

GNC.Control.LQR.Model = 'RLM';                  % LQR model
GNC.Control.LQR.Q = 1e2*eye(9);                 % Penalty on the state error
GNC.Control.LQR.M = eye(3);                     % Penalty on the control effort
GNC.Control.LQR.Reference = Sn(1,1:3);          % Reference operting point

GNC.Tmax = 1e-3 / (4*pi^2*Lem/Tc^2);            % Maximum available acceleration

% Augmented initial conditions 
int = zeros(1,3);                               % Integral of the relative position
s0 = [S1(end,:) int];                           % Initial conditions

dt = 1e-3;                                      % Continuous time step
tf(2) = 1.5*pi;                                 % Approach TOF
tspan = 0:dt:tf(2);                             % Integration time span

[~, Sg_lqr] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);

% Relative formation flying trajectory
tspan = 0:dt:pi/3;
tf(2) = tf(2) + pi/3;
Sg_for = repmat(Sg_lqr(end,n+1:2*n), length(tspan), 1);
s0 = [Sg_lqr(end,1:2*n) 0 0 0];

order = 10; 
[Cp, Cv, ~] = CTR_guidance(order, tspan, Sg_for);

GNC.Algorithms.Guidance = 'CTR';                             % Guidance algorithm
GNC.Guidance.CTR.Order = order;                              % Order of the approximation
GNC.Guidance.CTR.TOF = tspan(end);                           % Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	         % Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;                  % Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.IntegralCoefficients = zeros(size(Cp));     % Coefficients of the Chebyshev approximation

[~, Sg_lqr2] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);

Sg_lqr = [Sg_lqr(:,1:2*n); Sg_lqr2(:,1:2*n)];

% Final guidance regression 
tspan = 0:dt:tf(2);                                          % Complete phase time span 

order = 50; 
[Cp, Cv, Cg] = CTR_guidance(order, tspan, Sg_lqr(1:end-1,n+1:2*n));

GNC.Algorithms.Guidance = 'CTR';                             % Guidance algorithm
GNC.Guidance.CTR.Order = order;                              % Order of the approximation
GNC.Guidance.CTR.TOF = tspan(end);                           % Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	         % Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;                  % Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;              % Coefficients of the Chebyshev approximation

% Control 
GNC.Navigation.NoiseVariance = 1/Lem;

GNC.Algorithms.Control = 'SMC';                              % Control algorithm

GNC.Control.SMC.Parameters = [1.000000000000000 0.985999332287318 0.006010671478548 0.013227007322678]; 

s0 = S1(end,:);                                              % Phase initial conditions 

tic;
[~, S2] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
Time(2) = toc; 

% Control law
[~, ~, u_smc] = GNC_handler(GNC, S2(:,1:n), S2(:,7:end), tspan);

% Error in time 
[e_smc, merit_smc] = figures_merit(tspan, S2(:,7:12));

% Control integrals
effort_smc = control_effort(tspan, u_smc, true);

% Range reduction 
range(2) = 1-norm(S2(end,7:9))/norm(S2(1,7:9));

%% Final rendezvous: MPC-MISG
% Scheme setup
tol = [1e-7 1e-3];             % Convergence toleranes 
N = 100;                       % Number of impulses
method = 'MPC';                % Solver method       

tf(3) = 0.6;                   % Final phase TOF 

s0 = S2(end,:);                % Phase initial conditions

% Controller scheme
tic
[tspan_misg, S3, u_misg, state] = MISG_control(mu, tf(3), s0, method, N, tol);  
Time(3) = toc;

% Total maneuver metrics 
effort_misg = control_effort(tspan_misg, u_misg, true);

% Error in time 
[e_misg, merit_misg] = figures_merit(tspan_misg, S3(:,7:12));

%% Save results 
if (false)
    load JWST_mission;
else
    save JWST_mission;
end

%% Final results 
% Complete trajectory 
St = [S1(:,1:2*n); S2(2:end,1:2*n); S3(2:end,1:2*n)];

% Total integration time
tspan = [0:5e-2:tf(1)  (1e-3:1e-3:tf(2)-1e-3)+tf(1) tf(2)+tspan_misg];                                                    

%% Plotting
figure 
view(3) 
hold on 
t = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'k');
c = plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'k-*', 'MarkerIndices', 1:200:size(Sn,1));
plot3(St(:,1)+St(:,7), St(:,2)+St(:,8), St(:,3)+St(:,9), 'r', 'LineWidth', 0.5); 
scatter3(St(end,1)+St(end,7), St(end,2)+St(end,8), St(end,3)+St(end,9), 'r', 'filled');
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
hold off
text(L(1,Ln)+1e-4, L(2,Ln), 1e-6, '$L_2$');
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
legend('Target orbit', 'Loiter orbit', 'Chaser motion', 'Location', 'northwest');
grid on;

% Configuration space evolution
figure
subplot(1,2,1)
hold on
plot(tspan(1:size(St,1)), St(:,7)); 
plot(tspan(1:size(St,1)), St(:,8)); 
plot(tspan(1:size(St,1)), St(:,9)); 
hold off
xlabel('$t$');
ylabel('$\mathbf{\rho}$');
grid on;
legend('$x$', '$y$', '$z$');
subplot(1,2,2)
hold on
plot(tspan(1:size(St,1)), St(:,10)); 
plot(tspan(1:size(St,1)), St(:,11)); 
plot(tspan(1:size(St,1)), St(:,12)); 
hold off
xlabel('$t$');
ylabel('$\dot{\mathbf{\rho}}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');

% Configuration space evolution (Phase I)
figure
subplot(1,2,1)
hold on
plot(tspan(1:size(St,1)), St1(:,7)); 
plot(tspan(1:size(St,1)), St1(:,8)); 
plot(tspan(1:size(St,1)), St1(:,9)); 
hold off
xlabel('$t$');
ylabel('$\mathbf{\rho}$');
grid on;
legend('$x$', '$y$', '$z$');
subplot(1,2,2)
hold on
plot(tspan(1:size(St,1)), St1(:,10)); 
plot(tspan(1:size(St,1)), St1(:,11)); 
plot(tspan(1:size(St,1)), St1(:,12)); 
hold off
xlabel('$t$');
ylabel('$\dot{\mathbf{\rho}}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');

% Configuration space evolution (Phase II)
figure
subplot(1,2,1)
hold on
plot(tspan(1:size(St,1)), St2(:,7)); 
plot(tspan(1:size(St,1)), St2(:,8)); 
plot(tspan(1:size(St,1)), St2(:,9)); 
hold off
xlabel('$t$');
ylabel('$\mathbf{\rho}$');
grid on;
legend('$x$', '$y$', '$z$');
subplot(1,2,2)
hold on
plot(tspan(1:size(St,1)), St2(:,10)); 
plot(tspan(1:size(St,1)), St2(:,11)); 
plot(tspan(1:size(St,1)), St2(:,12)); 
hold off
xlabel('$t$');
ylabel('$\dot{\mathbf{\rho}}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');

% Configuration space evolution (Phase III)
figure
subplot(1,2,1)
hold on
plot(tspan_misg, St3(:,7)); 
plot(tspan_misg, St3(:,8)); 
plot(tspan_misg, St3(:,9)); 
hold off
xlabel('$t$');
ylabel('$\mathbf{\rho}$');
grid on;
legend('$x$', '$y$', '$z$');
subplot(1,2,2)
hold on
plot(tspan_misg, St3(:,10)); 
plot(tspan_misg, St3(:,11)); 
plot(tspan_misg, St3(:,12));
hold off
xlabel('$t$');
ylabel('$\dot{\mathbf{\rho}}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');


if (false)
    dh = 50; 
    steps = fix(size(St,1)/dh);
    M = cell(1,steps);
    h = figure;
    filename = 'webb.gif';
    view([40 20])
    hold on
    plot3(flip(St0(:,1)), flip(St0(:,2)), flip(St0(:,3)), '.r', 'Linewidth', 0.1);
    plot3(St(:,1), St(:,2), St(:,3), '.-b', 'Linewidth', 0.1);
    plot3(St4(:,1)+St4(:,7), St4(:,2)+St4(:,8), St4(:,3)+St4(:,9), '.r', 'Linewidth', 0.1); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
    scatter3(1-mu, 0, 0, 'k', 'filled');
    text(L(1,Ln)-2e-3, L(2,Ln), 0, '$L_2$');
    text(1-mu-1e-3, 0, 1e-3, '$M_2$');
    grid on;
    title('Rendezvous simulation');
    
    for i = 1:dh:size(St,1)
        T = scatter3(St(i,1), St(i,2), St(i,3), 20, 'b', 'filled'); 
        V = scatter3(St(i,1)+St(i,7), St(i,2)+St(i,8), St(i,3)+St(i,9), 20, 'r', 'filled');
        drawnow;
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); 
        if (i == 1) 
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1e-3); 
        else 
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1e-3); 
        end 
        delete(T); 
        delete(V);
    end
    hold off
end