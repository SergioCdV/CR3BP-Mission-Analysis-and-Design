%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 10/07/21 % 

%% GNC 11: Complete rendezvous mission example 4 %% 
% This script provides an interface to test the general control scheme for a rendezvous, docking and undocking mission. 

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L2 in the
% Sun-Earth system is analyzed.

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
mu = 3.003e-6;                      %Earth-Moon reduced gravitational parameter
Lem = 149597870700;                 %Mean distance from the Earth to the Moon

%% Initial conditions and halo orbit computation %%
%Initial conditions
L = libration_points(mu);                                   %System libration points
Az = 195e6;                                                 %Orbit amplitude out of the synodic plane. Play with it! 
Az = dimensionalizer(384400e3, 1, 1, Az, 'Position', 0);    %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param = [1 Az Ln gamma m];                                  %Halo orbit parameters (-1 being for southern halo)

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

%Generate the orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
Sn = target_orbit.Trajectory; 

%% Compute the tranfer trajectory
%Parking orbit definition 
R = [1-mu; 0; 0];                       %Primary location in the configuration space
branch = 'R';                           %Manifold branch to globalize
map = 'Secondary primary';              %Poincaré map to use
event = @(t,s)sp_crossing(t,s,mu);      %Integration event

hd = dimensionalizer(Lem, 1, 1, 2000e3, 'Position', 0);                  %Parking orbit altitude

%Integrate the stable manifold backwards and check if it intersects the whereabouts of the parking orbit
manifold = 'S';                                                          %Integrate the stable manifold
seed = Sn;                                                               %Periodic orbit seed
tspan = 0:1e-3:target_orbit.Period;                                      %Original integration time
rho = 20;                                                                %Density of fibres to analyze
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

%% First phase: long-range rendezvous using the MPC approach
%Time of flight 
tf(1) = 0.6;                                  %Nondimensional maneuver end time 

%Set up of the optimization
method = 'NPL';                               %Method to solve the problem
core = 'Linear';                              %Number of impulses
TOF = tf(1);                                  %Time of flight
cost_function = 'Position';                   %Target a position rendezvous

%Thruster characteristics 
Tmin = -0.1;                                  %Minimum thrust capability (in velocity impulse)
Tmax = 0.1;                                   %Maximum thrust capability (in velocity impulse)

%Main computation 
tic
[St1, dV, ~] = MPC_control(mu, cost_function, Tmin, Tmax, TOF, s0, core, method);
toc

%Control integrals
effort_mpc = control_effort(tspan, dV, true);

%% GNC: SDRE/LQR control law
%Phase definition 
tf(2) = pi; 
tspan = 0:dt:tf(2); 

%Guidance 
A = dimensionalizer(Lem, 1, 1, 100, 'Position', 0)*[1 1 0];
Sg = [A.*ones(length(tspan),3) zeros(length(tspan),3)];
order = 5; 
[Cp, Cv, Cg, Ci] = CTR_guidance(order, tspan, Sg);

model = 'RLM';
GNC.Algorithms.Guidance = 'CTR';                    %Guidance algorithm
GNC.Algorithms.Navigation = '';                     %Navigation algorithm
GNC.Algorithms.Control = 'SDRE';                    %Control algorithm
GNC.Guidance.Dimension = 9;                         %Dimension of the guidance law
GNC.Control.Dimension = 3;                          %Dimension of the control law

GNC.System.mu = mu;                                 %Systems's reduced gravitational parameter
GNC.System.Libration = [Ln gamma];                  %Libration point ID

GNC.Control.SDRE.Model = model;                     %SDRE model
GNC.Control.SDRE.Q = 2*eye(9);                      %Penalty on the state error
GNC.Control.SDRE.M = eye(3);                        %Penalty on the control effort

GNC.Guidance.CTR.Order = order;                     %Order of the approximation
GNC.Guidance.CTR.TOF = tf(2);                       %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	%Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;         %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;     %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.IntegralCoefficients = Ci;         %Coefficients of the Chebyshev approximation

%Initial conditions 
int = zeros(1,3);                                   %Integral of the relative position
slqr0 = [St1(end,1:12) int];                        %Initial conditions

%Compute the trajectory
tic
[~, St2] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, slqr0, options);
toc

%Control law
[~, ~, u] = GNC_handler(GNC, St2(:,1:6), St2(:,7:end), tspan);

%Control integrals
efforr_sdre = control_effort(tspan, u, false);

%% Third phase: rendezvous with MI
%Phase definition 
tf(3) = 0.8; 

%Differential corrector set up
tol = 1e-8;                                   %Differential corrector tolerance

%Select impulsive times 
times = [0 tf(3)*rand(1,5)];                  %Times to impulse the spacecraft

%Compute the control law
impulses.Number = length(times);              %Number of impulses
impulses.Weights = eye(impulses.Number*3);    %Weightening matrix
impulses.Times = times;                       %Impulses times

cost = 'Position';                            %Cost function to target

%Compute the guidance law
tic
[St3, dV3, state] = MISS_control(mu, tf(3), St2(end,1:12), tol, cost, impulses);
toc

%Performance indices
effort_mi = control_effort(tspan, dV3, true); %Control effort made

%% Third phase: escape
%Phase definition
tf(4) = 1.5*pi;                                 %End of the escape 

%Controller definition
constraint.Constrained = false;                 %No constraints on the maneuver
constraint.SafeDistance = 1e-5;                 %Safety distance at the collision time
constraint.Period = target_orbit.Period(end);   %Orbital period
constraint.Energy = false;                      %No energy constraint

tic
[St4, dV4, tm] = FMSC_control(mu, 0.1, St3(end,1:12), tol, constraint, 'Center');
toc

%Re-integrate trajectory
tspan = 0:dt:tf(4)-0.1;
[~, St4b] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, St4(end,1:12), options);
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
text(L(1,Ln)+1e-4, L(2,Ln), 1e-3, '$L_2$');
text(1-mu-5e-4, 0, 1e-3, '$M_2$');
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
indexOfInterest = (tspan < 0.98*sum(tf(3))) & (tspan > 0); 
hold on
plot(tspan(indexOfInterest), St(indexOfInterest, 7))  
plot(tspan(indexOfInterest), St(indexOfInterest, 8))  
plot(tspan(indexOfInterest), St(indexOfInterest, 9))  
hold off
axis tight

subplot(1,2,2)
axes('position', [0.62 .30 .25 .25])
box on
indexOfInterest = (tspan < 0.98*sum(tf(3))) & (tspan > 0); 
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