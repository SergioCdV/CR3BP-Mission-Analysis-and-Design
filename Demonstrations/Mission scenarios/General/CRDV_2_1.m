%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 10/07/21 % 

%% GNC 11: Complete rendezvous mission example 2 %% 
% This script provides an interface to test the general control scheme for a rendezvous, docking and undocking mission. 

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

%CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
Lem = 384400e3;                     %Mean distance from the Earth to the Moon

%% Initial conditions and halo orbit computation %%
%Initial conditions
L = libration_points(mu);                                   %System libration points
Az = 195e6;                                                 %Orbit amplitude out of the synodic plane. Play with it! 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param = [-1 Az Ln gamma m];                                 %Halo orbit parameters (-1 being for southern halo)

%Correction parameters 
dt = 1e-3;                                                  %Time step to integrate converged trajectories
maxIter = 20;                                               %Maximum allowed iterations in the differential correction schemes
tol = 1e-10;                                                %Differential correction tolerance 
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                    %Number of orbits to continuate
direction = -1;                                             %Direction to continuate (to the Earth)
   
%% Functions
%Compute the NRHO
[halo_seed, haloT] = object_seed(mu, param, 'Halo');        %Generate a halo orbit seed

%Continuation procedure 
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, haloT};                       %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
setup = [mu maxIter tol direction];                         %General setup

[target_orbit, state_energy] = continuation(num, method, algorithm, object, corrector, setup);

%Generate the NRHO
s0 = [target_orbit.Seeds(end,:).'; reshape(eye(6), [36 1])];
tspan = 0:dt:target_orbit.Period(end);
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);
Sn = S; 

%% Compute the tranfer trajectory
%Parking orbit definition 
R = [1-mu; 0; 0];                       %Primary location in the configuration space
branch = 'R';                           %Manifold branch to globalize
map = 'Secondary primary';              %Poincaré map to use
event = @(t,s)sp_crossing(t,s,mu);      %Integration event

hd = dimensionalizer(Lem, 1, 1, 2000e3, 'Position', 0);                  %Parking orbit altitude

%Integrate the stable manifold backwards and check if it intersects the whereabouts of the parking orbit
manifold = 'S';                                                          %Integrate the stable manifold
seed = S;                                                                %Periodic orbit seed
tspan = 0:1e-3:target_orbit.Period;                                      %Original integration time
rho = 50;                                                                %Density of fibres to analyze
S = invariant_manifold(mu, Ln, manifold, branch, seed, rho, tspan, map); %Initial trajectories

%Relative distance to the primary of interest
distance = zeros(rho,1);    
for i = 1:rho
    %Distance to the orbital altitude
    distance(i) = norm(shiftdim(S.Trajectory(i,S.ArcLength(i),1:3))-R)-hd;  
end

[~, index] = sort(distance);                            %Select the closest manifold to the parking orbit
s0 = shiftdim(S.Trajectory(index(1),1,:));              %Initial conditions to correct
Phi = eye(n);                                           %Initial STM 
Phi = reshape(Phi, [n^2 1]);                            %Initial STM 
sHalo = seed(S.Index(index(1)),1:n).';                  %Halo insertion point
s0 = [s0(1:3); s0(4:6); Phi];                           %Initial conditions
TOF = S.TOF(index(1));                                  %Time of flight
tspan = TOF:-dt:0;                                                                  %Integration time
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', event);             %Integration tolerances  
[~, St0] = ode113(@(t,s)cr3bp_equations(mu, 1, true, t, s), tspan, s0, options);    %Natural trajectory

%% Modelling in the synodic frame %%
r_t0 = Sn(mod(size(St0,1), size(Sn,1)), 1:6);               %Initial target conditions
r_c0 = St0(1,1:6);                                          %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Insertion plot 
figure(1)
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3));
plot3(St0(:,1), St0(:,2), St0(:,3));
scatter3(r_c0(1), r_c0(2), r_c0(3))
scatter3(r_t0(1), r_t0(2), r_t0(3))
hold off
legend('Target orbit', 'Target at HOI', 'Chaser at HOI'); 
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
title('Reconstruction of the natural chaser motion');

%% First phase: long-range rendezvous using the TITA approach
%Time of flight 
tf(1) = 0.6;                                %Nondimensional maneuver end time 
sd = zeros(1,3);                            %Desired state of the system

%Differential corrector set up
tol = 1e-5;                                 %Differential corrector tolerance

%Cost function matrices
penalties.R = eye(3);                       %Penalty on the impulse
penalties.Q = eye(6);                       %Penalty on the state error
penalties.M  = 0.1*eye(6);                  %Penalty on the state noise

%Select measuring times 
target_points.Noise = true;                 %Boolean to account for state noise
target_points.Times = tf(1)*rand(1,3);      %Times to measure the state noise

thruster_model.Sigma = 0.01;                %Velocity noise dependance on the velocity impulse
thruster_model.Rotation = eye(3);           %Rotational misalignment of the thrusters

%Cost function 
cost_function = 'Position';                 %Cost function to target
two_impulsive = true;                       %Two-impulsive rendezvous boolean

%TITA controller
tic
[St1, dV, state] = TITA_control(mu, tf(1), s0, tol, cost_function, sd, two_impulsive, penalties, target_points, thruster_model);
toc

%Control effort 
tspan = 0:dt:tf(1);
effort_tita = control_effort(tspan, dV, true);

%% Second phase: close-range rendezvous
%Phase definition 
tf(2) = pi;                                 %Timeline

%GNC algorithms definition 
GNC.Algorithms.Guidance = '';               %Guidance algorithm
GNC.Algorithms.Navigation = '';             %Navigation algorithm
GNC.Algorithms.Control = 'SMC';             %Control algorithm
GNC.Algorithms.Solver = 'Encke';            %Dynamics vector field to be solved
GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                  %Dimension of the control law
GNC.System.mu = mu;                         %System reduced gravitational parameter

%Controller parameters
%GNC.Control.SMC.Parameters = [1 SMC_optimization(mu, 'L2', St1(end,1:12), tf(2))]; 
GNC.Control.SMC.Parameters = [1 1.000000000000000 0.985999332287318 0.006010671478548 0.013227007322678]; 

%Integration time 
tspan = 0:dt:tf(2); 

%Initial conditions 
s0 = St1(end,:); 

%Re-integrate trajectory
tic
[~, St2] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Control effort
[~, ~, u] = GNC_handler(GNC, St2(:,1:6), St2(:,7:12), tspan); 
effort_smc1 = control_effort(tspan, u, false);

%% Third phase: formation-flying
%Phase definition 
tf(3) = pi; 
tspan = 0:dt:tf(3); 

%Guidance 
A = dimensionalizer(Lem, 1, 1, 1e6, 'Position', 0);
Sg = A*[cos(2*pi*tspan/tspan(end)).' sin(2*pi*tspan/tspan(end)).' zeros(length(tspan),4)];
order = 5; 
[Cp, Cv, Cg] = CTR_guidance(order, tspan, Sg);

%Reconstructed guidance trajectory
u = (2*tspan-(tspan(end)+tspan(1)))/(tspan(end)-tspan(1));          %Normalized time domain
T = CH_basis('first', order, u);                                    %Polynomial basis

%Error in the regression
p = Cp*T;                   %Position regression

%GNC algorithms definition 
GNC.Algorithms.Guidance = 'CTR';            %Guidance algorithm
GNC.Algorithms.Navigation = '';             %Navigation algorithm
GNC.Algorithms.Control = 'SMC';             %Control algorithm
FNC.Algorithms.Solver = 'Encke';            %Dynamics vector field to be solved
GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                  %Dimension of the control law
GNC.System.mu = mu;                         %System reduced gravitational parameter

%Controller parameters
%GNC.Control.SMC.Parameters = [1 SMC_optimization(mu, 'L2', St2(end,1:12), tf(2))]; 
GNC.Control.SMC.Parameters = [1 1.000000000000000 0.432562054680836 0.070603623964497 0.099843662546135]; 

%Guidance parameters 
GNC.Guidance.CTR.Order = order;                     %Order of the approximation
GNC.Guidance.CTR.TOF = tf(3);                       %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	%Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;         %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;     %Coefficients of the Chebyshev approximation

%Integration time 
tspan = 0:dt:tf(3); 

%Initial conditions 
s0 = St2(end,:); 

%Re-integrate trajectory
tic
[~, St3] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Control effort
[~, ~, u] = GNC_handler(GNC, St3(:,1:6), St3(:,7:12), tspan); 
effort_smc2 = control_effort(tspan, u, false);

%% Third phase: escape
%Phase definition
tf(4) = 2*pi;                                 %End of the escape 

%Controller definition
constraint.Constrained = false;               %No constraints on the maneuver
constraint.SafeDistance = 1e-5;               %Safety distance at the collision time
constraint.Period = target_orbit.Period(end); %Period of the target period
constraint.Energy = false;                    %No energy constraint

tic
[St4, dV4, tm] = FMSC_control(mu, 0.1, St3(end,1:12), tol, constraint, 'Center');
toc

%Re-integrate trajectory
tspan = 0:dt:tf(4)-0.1;
[~, St4b] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, St4(end,:), options);
St4 = [St4(1:end-1,:); St4b];

%Control effort 
effort_fmsc = control_effort(tspan, dV4, true);

%% Final results 
%Complete trajectory 
St = [St1(1:end-1,1:12); St2(1:end-1,1:12); St3(1:end-1,1:12); St4(1:end,1:12)];

%Total integration time
tspan = 0:dt:sum(tf);                                                    

%% Plotting
figure(3) 
view(3) 
hold on 
t = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'b', 'Linewidth', 0.1);
plot3(flip(St0(:,1)), flip(St0(:,2)), flip(St0(:,3)), 'r', 'Linewidth', 0.1);
plot3(St(:,1)+St(:,7), St(:,2)+St(:,8), St(:,3)+St(:,9), 'r', 'Linewidth', 0.1); 
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
scatter3(1-mu, 0, 0, 'k', 'filled');
hold off
text(L(1,Ln)-5e-3, L(2,Ln)+1e-3, 5e-3, '$L_2$');
text(1-mu+1e-3, 0, 5e-3, '$M_2$');
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
legend('Target orbit', 'Chaser trajectory', 'Location', 'northeast');
grid on;
title('Mission trajectory in the absolute configuration space');

%Configuration space evolution
figure(4)
subplot(1,2,1)
hold on
plot(tspan(1:size(St,1)), St(:,7)); 
plot(tspan(1:size(St,1)), St(:,8)); 
plot(tspan(1:size(St,1)), St(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('$x$', '$y$', '$z$');
title('Relative position evolution');

subplot(1,2,2)
hold on
plot(tspan(1:size(St,1)), St(:,10)); 
plot(tspan(1:size(St,1)), St(:,11)); 
plot(tspan(1:size(St,1)), St(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinate');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Relative velocity evolution');

subplot(1,2,1)
axes('position', [.17 .52 .25 .25])
box on
indexOfInterest = (tspan < 0.98*sum(tf(1:4))) & (tspan > sum(tf(1))); 
hold on
plot(tspan(indexOfInterest), St(indexOfInterest, 7))  
plot(tspan(indexOfInterest), St(indexOfInterest, 8))  
plot(tspan(indexOfInterest), St(indexOfInterest, 9))  
hold off
axis tight

subplot(1,2,2)
axes('position', [0.62 .30 .25 .25])
box on
indexOfInterest = (tspan < 0.98*sum(tf(1:4))) & (tspan > sum(tf(1))); 
hold on
plot(tspan(indexOfInterest), St(indexOfInterest, 10))  
plot(tspan(indexOfInterest), St(indexOfInterest, 11))  
plot(tspan(indexOfInterest), St(indexOfInterest, 12))  
hold off
axis tight

if (false)
    dh = 50; 
    steps = fix(size(St,1)/dh);
    M = cell(1,steps);
    h = figure;
    filename = 'nhro.gif';
    view([37 20])
    hold on
    plot3(flip(St0(:,1)), flip(St0(:,2)), flip(St0(:,3)), '.b', 'Linewidth', 0.1);
    plot3(St(:,1), St(:,2), St(:,3), '.-b', 'Linewidth', 0.1);
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
    scatter3(1-mu, 0, 0, 'k', 'filled');
    text(L(1,Ln)-2e-2, L(2,Ln), 0, '$L_2$');
    text(1-mu+1e-3, 0, 5e-3, '$M_2$');
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