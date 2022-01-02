%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 29/12/21 % 

%% JGCD Gateway Mission %% 
% This script provides an interface to test the general control scheme for
% a rendezvous, docking and undocking mission in a re-supply mission for
% the Gateway 

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
R = [-mu; 0; 0];                        %Primary location in the configuration space
branch = 'R';                           %Manifold branch to globalize
map = 'First primary';                  %Poincaré map to use
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
penalties.Q = 2*eye(6);                     %Penalty on the state error
penalties.M  = 0.1*eye(6);                  %Penalty on the state noise

%Select measuring times 
target_points.Noise = true;                 %Boolean to account for state noise

thruster_model.Sigma = 0.01;                %Velocity noise dependance on the velocity impulse
thruster_model.Rotation = eye(3);           %Rotational misalignment of the thrusters

%Cost function 
cost_function = 'Position';                 %Cost function to target
two_impulsive = true;                       %Two-impulsive rendezvous boolean

%TITA controller first stage
target_points.Times = 0.8*tf(1)*rand(1,3);  %Times to measure the state noise
tic
[Staux1, dVaux1, ~] = TITA_control(mu, 0.8*tf(1), s0, cost_function, zeros(1,3), two_impulsive, penalties, target_points, thruster_model, tol);
toc

%Multi-impulse stage
tol = 1e-8;                                   %Differential corrector tolerance

%Select impulsive times 
times = [0.2*tf(1)*rand(1,5)];                %Times to impulse the spacecraft

%Compute the control law
impulses.Number = length(times);              %Number of impulses
impulses.Weights = eye(impulses.Number*3);    %Weightening matrix
impulses.Times = times;                       %Impulses times

cost = 'Position';                            %Cost function to target

%Controller scheme
tic
[Staux2, dVaux2, state] = MISS_control(mu, 0.2*tf(1), Staux1(end,1:2*n), tol, cost, impulses);
toc
toc

St1 = [Staux1(1:end,1:2*n); Staux2(2:end,1:2*n)];
dV = [dVaux1 dVaux2]; 

%Control effort 
tspan = 0:dt:tf(1);
effort_first = control_effort(tspan, dV, true);

%% Second phase: close-range rendezvous
%Phase definition 
tf(2) = pi/14;                              %Timeline

%GNC algorithms definition 
GNC.Algorithms.Guidance = '';               %Guidance algorithm
GNC.Algorithms.Navigation = '';             %Navigation algorithm
GNC.Algorithms.Control = 'SMC';             %Control algorithm
GNC.Algorithms.Solver = 'Encke';            %Dynamics vector field to be solved
GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                  %Dimension of the control law
GNC.System.mu = mu;                         %System reduced gravitational parameter

GNC.Navigation.NoiseVariance = dimensionalizer(Lem, 1, 1, 50, 'Position', 0);

%Controller parameters
GNC.Control.SMC.Parameters = [1 1.000000000000000 0.985999332287318 0.006010671478548 0.013227007322678]; 

%Integration time 
tspan = 0:dt:tf(2); 

%Initial conditions 
s0 = St1(end,:); 

%Re-integrate trajectory
tic
[~, St2] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Control effort
[~, ~, u] = GNC_handler(GNC, St2(:,1:6), St2(:,7:12), tspan); 
effort_second = control_effort(tspan, u, false);

%% Third phase: formation-flying
% Guidance (Lissajous trajectory around the target)
%Phase definition 
tf(3) = pi/2; 
tspan = 0:dt:tf(3); 

%Target trajectory
s0 = St2(end,1:6);
[~, Staux] = ode113(@(t,s)cr3bp_equations(mu, 1, false, t, s), tspan, s0, options);   

%Compute the Frenet-Serret
T = zeros(size(Staux,1), 3, 3);    %Preallocation of the frame
dT = T;                            %Preallocation of the angular velocity of the frame
for i = 1:size(Staux,1)
    T(i,:,:) = frenet_triad(mu, Staux(i,1:6));                 %Frenet-Serret triad
end

for i = 1:size(Staux,1)
    if (i == 1)
        dT(i,:,:) = (T(i+1,:,:)-T(size(Staux,1),:,:))/(2*dt);  %Frenet-Serret triad derivative
    elseif (i == size(Staux,1))
        dT(i,:,:) = (T(1,:,:)-T(i,:,:))/(2*dt);                %Frenet-Serret triad derivative
    else
        dT(i,:,:) = (T(i+1,:,:)-T(i,:,:))/(2*dt);              %Frenet-Serret triad derivative
    end
end

%Generate an 1 km helix around in the Frenet-Serret XY plane for the given
A = dimensionalizer(Lem, 1, 1, 50, 'Position', 0);
omega = 4*pi/tspan(end);
Sg_frenet(:,1:3) = [dt*tspan.'/100 A*cos(omega*tspan).' A*sin(omega*tspan).'];
Sg_frenet(:,4:6) = zeros(length(tspan),3);
Sg = Sg_frenet; 
for i = 1:size(Sg,1)
    omega = zeros(1,3);
    for j = 1:size(T,3)
        omega = omega + (1/2)*cross(T(i,:,j),dT(i,:,j));     %Darboux vector
    end
    Sg(i,1:3) = (shiftdim(T(i,:,:)).'*Sg(i,1:3).').';        %Synodic position
    Sg(i,4:6) = (shiftdim(T(i,:,:)).'*Sg(i,4:6).').';        %Synodic velocity
    Sg(i,4:6) = Sg(i,4:6) + cross(omega,Sg(i,4:6));          %Synodic velocity Coriolis correction
end

%Compute the trajectory as a Chebyshev analytical expression
order = 100; 
[Cp, Cv, Cg] = CTR_guidance(order, tspan, Sg);

%Reconstructed guidance trajectory
T = zeros(order, length(tspan));                                    %Preallocation of the polynomial basis
u = (2*tspan-(tspan(end)+tspan(1)))/(tspan(end)-tspan(1));          %Normalized time domain

for i = 1:length(tspan)
    T(:,i) = chebyshev('first', order, u(i));
end

%Error in the regression
p = Cp*T;                   %Position regression
v = Cv*T;                   %Velocity regression
Sr = Staux+[p.' v.'];       %Regress the phase space trajectory

%GNC algorithms definition 
GNC.Algorithms.Guidance = 'CTR';            %Guidance algorithm
GNC.Algorithms.Navigation = '';             %Navigation algorithm
GNC.Algorithms.Control = 'SMC';             %Control algorithm
FNC.Algorithms.Solver = 'Encke';            %Dynamics vector field to be solved
GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                  %Dimension of the control law
GNC.System.mu = mu;                         %System reduced gravitational parameter

GNC.Navigation.NoiseVariance = dimensionalizer(Lem, 1, 1, 1, 'Position', 0);

%Controller parameters
%GNC.Control.SMC.Parameters = [1 SMC_optimization(mu, 'L2', St2(end,1:12), tf(2))]; 
GNC.Control.SMC.Parameters = [1 1.000000000000000 0.432562054680836 0.070603623964497 0.099843662546135]; 

%Guidance parameters 
GNC.Guidance.CTR.Order = order;                     %Order of the approximation
GNC.Guidance.CTR.TOF = tf(3);                       %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	%Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;         %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;     %Coefficients of the Chebyshev approximation

%Initial conditions 
s0 = St2(end,:); 

%Re-integrate trajectory
tic
[~, St3] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Control effort
[~, ~, u] = GNC_handler(GNC, St3(:,1:6), St3(:,7:12), tspan); 
effort_third = control_effort(tspan, u, false);

%% Fourth phase: formation-flying
%Integration time 
tf(4) = pi/4;

%Impulsive first approach
tic
[St4aux, dV4, state] = TISS_control(mu, 0.1*tf(4), St3(end,1:2*n), 1e-10, 'Position', true);  %Controller scheme
toc

effort_fourth1 = control_effort(NaN, dV4, true);

%GNC algorithms definition 
GNC.Algorithms.Guidance = '';               %Guidance algorithm
GNC.Algorithms.Navigation = '';             %Navigation algorithm
GNC.Algorithms.Control = 'SDRE';            %Control algorithm
FNC.Algorithms.Solver = 'Encke';            %Dynamics vector field to be solved
GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                  %Dimension of the control law
GNC.System.mu = mu;                         %System reduced gravitational parameter

GNC.Navigation.NoiseVariance = dimensionalizer(Lem, 1, 1, 0, 'Position', 0);

%Controller paramters
GNC.System.Libration = [Ln gamma];          %Libration point ID
GNC.Control.SDRE.Model = 'RLM';             %SDRE model
GNC.Control.SDRE.Q = 2*eye(9);              %Penalty on the state error
GNC.Control.SDRE.M = eye(3);                %Penalty on the control effort

%Initial conditions 
tspan = 0:dt:0.9*tf(4); 
s0 = [St4aux(end,1:2*n) zeros(1,3)]; 

%Re-integrate trajectory
tic
[~, St4] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Compute the trajectory as a Chebyshev analytical expression
order = 50; 
[Cp, Cv, Cg] = CTR_guidance(order, tspan, St4(:,n+1:2*n));

%GNC algorithms definition 
GNC.Algorithms.Guidance = 'CTR';            %Guidance algorithm
GNC.Algorithms.Navigation = '';             %Navigation algorithm
GNC.Algorithms.Control = 'SMC';             %Control algorithm
FNC.Algorithms.Solver = 'Encke';            %Dynamics vector field to be solved
GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                  %Dimension of the control law
GNC.System.mu = mu;                         %System reduced gravitational parameter

GNC.Navigation.NoiseVariance = dimensionalizer(Lem, 1, 1, 10, 'Position', 0);

%Controller parameters
GNC.Control.SMC.Parameters = [1 1.000000000000000 0.432562054680836 0.070603623964497 0.099843662546135]; 

%Guidance parameters 
GNC.Guidance.CTR.Order = order;                     %Order of the approximation
GNC.Guidance.CTR.TOF = tf(4);                       %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	%Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;         %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;     %Coefficients of the Chebyshev approximation 

%Re-integrate trajectory
tic
[~, St4] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0(1:2*n), options);
toc

%Control effort
[~, ~, u] = GNC_handler(GNC, St4(:,1:6), St4(:,7:12), tspan); 
effort_fourth = control_effort(tspan, u, false);

%Final trajectory 
St4 = [St4aux(:,1:2*n); St4];

%% Final departure
%Phase definition
tf(5) = 0.1;                             %End of the escape 
tspan = 0:dt:tf(5);

%Controller definition
constraint.Constrained = false;                 %No constraints on the maneuver
constraint.SafeDistance = 1e-5;                 %Safety distance at the collision time
constraint.Period = target_orbit.Period(end);   %Orbital Period
constraint.Energy = true;                       %Energy constraint

tic
[St5, dV5, tm] = FMSC_control(mu, tf(5), St4(end,1:12), 1e-8, constraint, 'Center');
toc 

%Control effort 
effort_fifth = control_effort(tspan, dV5, true);

tf(6) = 2;
tspan = 0:dt:tf(6);
[~, St6] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, St5(end,:), options);

%% Final results 
%Complete trajectory 
St = [St1(1:end-1,1:2*n); St2(2:end,1:2*n); St3(2:end,1:2*n); St4(2:end,1:2*n); St5(2:end,1:2*n); St6(2:end,1:2*n)];

%Total integration time
tspan = 0:dt:sum(tf);                                                    

%% Plotting
figure(3) 
view(3) 
hold on 
t = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'b*', 'MarkerIndices', 1:200:size(Sn,1));
plot3(flip(St0(:,1)), flip(St0(:,2)), flip(St0(:,3)), 'r');
plot3(St(:,1)+St(:,7), St(:,2)+St(:,8), St(:,3)+St(:,9), 'r'); 
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
axes('position', [.22 .28 .2 .2])
box on
indexOfInterest = (tspan < 0.98*sum(tf(1:5))) & (tspan > sum(tf(2))); 
hold on
plot(tspan(indexOfInterest), St(indexOfInterest, 7))  
plot(tspan(indexOfInterest), St(indexOfInterest, 8))  
plot(tspan(indexOfInterest), St(indexOfInterest, 9))  
hold off
axis tight
ax = gca; 
ax.YAxis.Exponent = 0;

subplot(1,2,2)
axes('position', [0.67 .50 .18 .2])
box on
indexOfInterest = (tspan < 0.98*sum(tf(1:5))) & (tspan > sum(tf(1))); 
hold on
plot(tspan(indexOfInterest), St(indexOfInterest, 10))  
plot(tspan(indexOfInterest), St(indexOfInterest, 11))  
plot(tspan(indexOfInterest), St(indexOfInterest, 12))  
hold off
axis tight
ax = gca; 
ax.YAxis.Exponent = 0;

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