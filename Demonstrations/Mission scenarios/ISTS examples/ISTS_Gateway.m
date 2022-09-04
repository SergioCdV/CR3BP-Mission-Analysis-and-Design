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

%% Insertion into the loitering orbit
%Parking orbit definition 
R = [-mu; 0; 0];                        %Primary location in the configuration space
branch = 'L';                           %Manifold branch to globalize
map = 'First primary';                  %Poincaré map to use
event = @(t,s)sp_crossing(t,s,mu);      %Integration event

hd = dimensionalizer(Lem, 1, 1, 35e6, 'Position', 0);                    %Parking orbit altitude (GEO belt)

%Integrate the stable manifold backwards and check if it intersects the whereabouts of the parking orbit
manifold = 'S';                                                          %Integrate the stable manifold
seed = loiter_orbit.Trajectory;                                          %Periodic orbit seed
tspan = 0:dt:loiter_orbit.Period;                                        %Original integration time
rho = 70;                                                                %Density of fibres to analyze
S = invariant_manifold(mu, Ln, manifold, branch, seed, rho, tspan, map); %Initial trajectories

manifold = S;

%Relative distance to the primary of interest
distance = zeros(rho,1);    
for i = 1:rho
    %Distance to the orbital altitude
    distance(i) = norm(shiftdim(S.Trajectory(i,S.ArcLength(i),1:3))-R)-hd;  
end

[~, index] = sort(distance);                                                        %Select the closest manifold to the parking orbit
s0 = shiftdim(S.Trajectory(index(1),1,:));                                          %Initial conditions to correct
sHalo = seed(S.Index(index(1)),1:n).';                                              %Halo insertion point
TOF = S.TOF(index(1));                                                              %Time of flight
tspan = TOF:-dt:0;                                                                  %Integration time
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', event);             %Integration tolerances  
[~, St0] = ode113(@(t,s)cr3bp_equations(mu, 1, false, t, s), tspan, s0, options);   %Natural insertion trajectory
tf(1) = TOF;

%% Resident stationkeeping
%Integration of the model
s0 = loiter_orbit.Trajectory(1,1:n);        %Initial resident conditions
s0 = [s0 reshape(eye(n), [1 n^2])];         %Initial resident conditions
tspan = 0:dt:loiter_orbit.Period;           %Integration time span 

[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options); 

%Guidance law coefficients
Sg = jacobi_constant(mu,S(1,1:n).');

GNC.Algorithms.Guidance = '';                       %Guidance algorithm
GNC.Algorithms.Navigation = '';                     %Navigation algorithm
GNC.Algorithms.Control = 'MLQR';                    %Control algorithm

GNC.Guidance.Dimension = 9;                         %Dimension of the guidance law
GNC.Control.Dimension = 3;                          %Dimension of the control law

GNC.System.mu = mu;                                 %Systems's reduced gravitational parameter
GNC.Control.MLQR.Q = eye(2);                        %Penalty on the state error
GNC.Control.MLQR.M = eye(3);                        %Penalty on the control effort
GNC.Control.MLQR.Reference = Sg;                    %Penalty on the control effort
GNC.Control.MLQR.Period = loiter_orbit.Period;      %Penalty on the control effort

%Floquet mode matrix
[E, lambda] = eig(reshape(Sn(end,n+1:end), [n n]));
GNC.Control.MLQR.FloquetModes = diag(log(diag(lambda))/loiter_orbit.Period(end));
for i = 1:size(E,2)
    E(:,i) = E(:,i)/lambda(i,i);
end
GNC.Control.MLQR.FloquetDirections = E; 

%Initial conditions 
k = dimensionalizer(Lem, 1, 1, 100, 'Position', 0);     %Noise gain
r_t0 = S(1,1:n);                                        %Initial guidance target conditions
s0 = r_t0-2*k*ones(1,n)+k*rand(1,n);                    %Noisy initial conditions

GNC.Navigation.NoiseVariance = 0;                       %Noise variance

%Compute the trajectory
s0 = [r_t0 s0-r_t0 reshape(eye(n), [1 n^2])];
tic
[~, Staux1] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan(tspan < 0.1), s0, options);
[~, Staux2] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan(size(Staux1,1):end), Staux1(end,:), options);
toc
Sresident = [Staux1; Staux2(2:end,:)];

%Control law
[~, ~, u] = GNCc_handler(GNC, Staux1(:,1:n), Staux1(:,n+1:end), tspan(tspan < 0.1));

%Control integrals
energy_sk = control_effort(tspan(tspan < 0.1), u, false);

total_cost(1) = 1000*1.025*norm(energy_sk(:,3));

%% Rendezvous with the resident
%Long-range rendezvous with MPC 
r_t0 = Sresident(500,1:n)+Sresident(500,n+1:2*n);   %Initial resident conditions                  
r_c0 = sHalo.';                                     %Initial transfer vehicle conditions 
rho0 = r_c0-r_t0;                                   %Initial relative conditions
s0 = [r_t0 rho0].';                                 %Initial conditions of the target and the relative state

%Set up of the optimization
method = 'NPL';                               %Method to solve the problem
core = 'Linear';                              %Number of impulses
TOF = 0.5;                                    %Time of flight
cost_function = 'Position';                   %Target a position rendezvous

tf(2) = TOF; 

%Thruster characteristics 
Tmin = -1e-1;                                 %Minimum thrust capability (in velocity impulse)
Tmax = 1e-1;                                  %Maximum thrust capability (in velocity impulse)

%Main computation 
tspan = 0:dt:TOF; 
tic
[St_mpc1, dV1, ~] = MPC_control(mu, cost_function, Tmin, Tmax, 0.8*TOF, s0, core, method);
toc

tic
[St_mpc2, dV2, ~] = MPC_control(mu, cost_function, Tmin, Tmax, 0.15*TOF, St_mpc1(end,1:2*n), core, method);
toc

tic
[St_mpc3, dV3, ~] = MPC_control(mu, cost_function, Tmin, Tmax, 0.03*TOF, St_mpc2(end,1:2*n), core, method);
toc

tic
[St_mpc4, dV4, ~] = MPC_control(mu, cost_function, Tmin, Tmax, 0.02*TOF, St_mpc3(end,1:2*n), core, method);
toc

St_mpc = [St_mpc1(:,1:2*n); St_mpc2(2:end,1:2*n); St_mpc3(2:end,1:2*n); St_mpc4(2:end,1:2*n)];
dV = [dV1 dV2 dV3 dV4];

%Control integrals
effort_mpc = control_effort(tspan, dV, true);

%Error in time 
[e_mpc, merit_mpc] = figures_merit(tspan, St_mpc);

%Total cost of the phase 
total_cost(2) = 1000*1.05*norm(effort_mpc(:,3));

%% Departure from the loitering orbit
%Parking orbit definition 
R = [-mu; 0; 0];                        %Primary location in the configuration space
branch = 'L';                           %Manifold branch to globalize
map = 'First primary';                  %Poincaré map to use
event = @(t,s)sp_crossing(t,s,mu);      %Integration event

hd = dimensionalizer(Lem, 1, 1, 35e6, 'Position', 0);                    %Parking orbit altitude (GEO belt)

%Integrate the stable manifold backwards and check if it intersects the whereabouts of the parking orbit
manifold = 'U';                                                          %Integrate the stable manifold
seed = loiter_orbit.Trajectory;                                          %Periodic orbit seed
tspan = 0:dt:loiter_orbit.Period;                                        %Original integration time
rho = 70;                                                                %Density of fibres to analyze
S = invariant_manifold(mu, Ln, manifold, branch, seed, rho, tspan, map); %Initial trajectories

%Relative distance to the primary of interest
distance = zeros(rho,1);    
for i = 1:rho
    %Distance to the orbital altitude
    distance(i) = norm(shiftdim(S.Trajectory(i,S.ArcLength(i),1:3))-R)-hd;  
end

[~, index] = sort(distance);                                                        %Select the closest manifold to the parking orbit
s0 = shiftdim(S.Trajectory(index(1),1,:));                                          %Initial conditions to correct
sHalo = seed(S.Index(index(1)),1:n).';                                              %Halo insertion point
TOF = S.TOF(index(1));                                                              %Time of flight
tspan = TOF:-dt:0;                                                                  %Integration time
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', event);             %Integration tolerances  
[~, Stf] = ode113(@(t,s)cr3bp_equations(mu, 1, false, t, s), tspan, s0, options);   %Natural insertion trajectory

tf(3) = TOF;

%% Phase 1 trajectory 
Sc1 = St_mpc;                                 %Transfer vehicle trajectory
tspan1 = 0:dt:sum(tf(2));                     %Time span

%% Gateway stationkeeping orbit 
%Integration of the model
s0 = [target_orbit.Seeds(end,:) reshape(eye(n), [1 n^2])];
tspan = 0:dt:target_orbit.Period(end);
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);
Sn = repmat(S, floor(2*pi/target_orbit.Period(end))+1, 1);     

%Guidance law coefficients
Sg = jacobi_constant(mu,Sn(1,1:n).');

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

%Initial conditions 
k = dimensionalizer(Lem, 1, 1, 100, 'Position', 0);     %Noise gain
r_t0 = Sn(1,1:n);                                       %Initial guidance target conditions
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

total_cost(3) = 1000^1.025*norm(energy_sk(:,3));

%% Heteroclinic insertion connection between the loitering and target orbit
%Manifold definition
TOF = pi;                       %Time of flight
rho = 90;                       %Manifold fibers to compute 

tf(4) = TOF; 

%Computation flags
position_fixed = false;         %Flag to determine a final target state
graphics = false;               %Flag to plot the manifolds

%Final target state
target_orbit.Trajectory = Sn;
target_orbit.TargetState = shiftdim(Sn(1,1:n)); 

%Connection itenerary
sequence = [1 2];               %Connection itenerary
branch = ['L' 'L'];             %Manifold branches to be propagated

%Trajectory design core
[Shr, ~] = HTRC_guidance(mu, sequence, branch, rho, target_orbit, loiter_orbit, TOF, position_fixed, graphics);

%Transfer trajectory
St_transfer = Shr.Trajectory; 

%% Rendezvous with the Gateway 
%Initial conditions
tf(5) = 0.9*pi/2;                      %Interception time
s0 = [Sgate(1,1:n)+Sgate(1,n+1:2*n) St_transfer(end,1:n)-(Sgate(1,1:n)+Sgate(1,n+1:2*n))];

%Initial two impulsive 
tic
[St4aux, dVti, ~] = TISS_control(mu, 0.1*tf(4), s0, 1e-10, 'Position', true);  %Controller scheme
toc
effort_ti = control_effort(NaN, dVti(1:end-1), true);

%SMC tracking 
GNC.Algorithms.Guidance = '';               	%Guidance algorithm
GNC.Algorithms.Navigation = '';                 %Navigation algorithm
GNC.Algorithms.Control = 'SMC';                 %Control algorithm
GNC.Guidance.Dimension = 9;                     %Dimension of the guidance law
GNC.Control.Dimension = 3;                      %Dimension of the control law
GNC.System.mu = mu;                             %System reduced gravitational parameters

GNC.Control.SMC.Parameters = [1.000000000000000 0.432562054680836 0.070603623964497 0.099843662546135];

%Re-integrate the trajectory
tspan = 0:dt:tf(5);
tic 
[~, S2_smc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, St4aux(1,1:2*n), options);
toc

%Error in time 
[e_smc, merit_smc] = figures_merit(tspan, S2_smc);

%Control law
[~, ~, u] = GNC_handler(GNC, S2_smc(:,1:6), S2_smc(:,7:12), tspan);    

%Control integrals
effort_smc = control_effort(tspan, u, false);

%Final impulsive approach using the multi-impulsive scheme
tol = 1e-8;                                   %Differential corrector tolerance

%Select impulsive times 
tf(6) = 0.8;                                  %Time to dock                
times = 0.9*tf(6)*rand(1,5);                  %Times to impulse the spacecraft

%Compute the control law
impulses.Number = length(times);              %Number of impulses
impulses.Weights = eye(impulses.Number*3);    %Weightening matrix
impulses.Times = times;                       %Impulses times

cost = 'Position';                            %Cost function to target

%Controller scheme
tic
[S2_miss1, dV1, ~] = MISS_control(mu, 0.9*tf(6), S2_smc(end,1:2*n), tol, cost, impulses);
toc

tic
times = 0.05*tf(6)*rand(1,5);                 %Times to impulse the spacecraft
impulses.Times = times;                       %Impulses times
[S2_miss2, dV2, state] = MISS_control(mu, 0.05*tf(6), S2_miss1(end,1:2*n), tol, cost, impulses);
toc

S2_miss = [S2_miss1(1:end,1:2*n); S2_miss2(2:end,1:2*n)];
dV = [dV1 dV2]; 

%Control effort 
effort_miss = control_effort(NaN, dV, true);

tic
[St4aux, dVti, ~] = TISS_control(mu, 0.05*tf(6), S2_miss(end,1:2*n), 1e-10, 'Position', true);  %Controller scheme
toc
effort_ti2 = control_effort(NaN, dVti, true);

%Totalc cost of the maneuver 
total_cost(5) = 1000*1.025*(norm(effort_smc(:,3))+norm(effort_ti(:,3))+norm(effort_ti2(:,3))+norm(effort_miss(:,3)));

%Complete rendezvous trajectory 
S2 = [S2_smc(:,1:2*n); S2_miss(2:end,1:2*n); St4aux(2:end,1:2*n)];

%% Heteroclinic departure connection between the loitering and target orbit
%Manifold definition
TOF = pi;                       %Time of flight
rho = 90;                       %Manifold fibers to compute 

tf(7) = TOF; 

%Computation flags
position_fixed = false;        %Flag to determine a final target state
graphics = true;              %Flag to plot the manifolds

%Final target state
target_orbit.Trajectory = Sn;
target_orbit.TargetState = shiftdim(S2(end,1:n)); 

%Connection itenerary
sequence = [2 1];               %Connection itenerary
branch = ['R' 'R'];             %Manifold branches to be propagated

%Trajectory design core
[Shr, dVhr] = HTRC_guidance(mu, sequence, branch, rho, target_orbit, loiter_orbit, TOF, position_fixed, graphics);
St_departure = Shr.Trajectory;

%Control integrals
energy_hr = control_effort(NaN, dVhr, true);

%% Phase 2 trajectory
Sc2 = S2(2:end,1:n)+S2(2:end,n+1:2*n);
tspan2 = 0:dt:sum(tf(5:6));

%% Results %% 
%Configuration space evolution
figure
subplot(1,2,1)
hold on
plot(tspan1, Sc1(:,7)); 
plot(tspan1, Sc1(:,8)); 
plot(tspan1, Sc1(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinates');
grid on;
legend('$x$', '$y$', '$z$');
title('Relative position in time');
ax = gca; 
ax.YAxis.Exponent = 0;
subplot(1,2,2)
hold on
plot(tspan1, Sc1(:,10)); 
plot(tspan1, Sc1(:,11)); 
plot(tspan1, Sc1(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinates');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Relative velocity in time');
ax = gca; 
ax.YAxis.Exponent = 0;

%Configuration space evolution
figure
subplot(1,2,1)
hold on
plot(tspan2, S2(:,7)); 
plot(tspan2, S2(:,8)); 
plot(tspan2, S2(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinates');
grid on;
legend('$x$', '$y$', '$z$');
title('Relative position in time');
ax = gca; 
ax.YAxis.Exponent = 0;
subplot(1,2,2)
hold on
plot(tspan2, S2(:,10)); 
plot(tspan2, S2(:,11)); 
plot(tspan2, S2(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinates');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Relative velocity in time');
ax = gca; 
ax.YAxis.Exponent = 0;

figure
view(3) 
hold on
r = plot3(Sresident(:,7)+Sresident(:,1), Sresident(:,8)+Sresident(:,2), Sresident(:,9)+Sresident(:,3), 'b', 'Marker', '+', 'MarkerIndices', 1:400:size(Sresident,1)); 
t = plot3(Sgate(:,7)+Sgate(:,1), Sgate(:,8)+Sgate(:,2), Sgate(:,9)+Sgate(:,3), 'k', 'Marker', '*', 'MarkerIndices', 1:400:size(Sresident,1));
g = plot3(Sc1(:,1)+Sc1(:,7), Sc1(:,2)+Sc1(:,8), Sc1(:,3)+Sc1(:,9), 'r');
h = plot3(Sc2(:,1), Sc2(:,2), Sc2(:,3), 'r');
plot3(St_transfer(:,1), St_transfer(:,2), St_transfer(:,3), 'r');
plot3(St_departure(:,1), St_departure(:,2), St_departure(:,3), 'r');
plot3(Stf(1:5000,1), Stf(1:5000,2), Stf(1:5000,3), 'r');
plot3(St0(1:2000,1), St0(1:2000,2), St0(1:2000,3), 'r');
scatter3(L(1,1), L(2,1), 0, 'k', 'filled');
scatter3(L(1,2), L(2,2), 0, 'k', 'filled');
scatter3(1-mu, 0, 0, 'k', 'filled');
hold off
text(L(1,1)+1e-3, L(2,1), 5e-3, '$L_1$');
text(L(1,2)+1e-3, L(2,2), 5e-3, '$L_2$');
text(1-mu+1e-3, 0, 5e-3, '$M_2$');
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Loiter orbit', 'Gateway orbit', 'Chaser trajectory', 'Location', 'northeast');
title('Mission trajectory in the absolute configuration space');

%Rendezvous animation 
if (true)
    dh = 100;
    W = figure(1);
    filename = 'gateway.gif';
    view([40 20]) 
    hold on
    r = plot3(Sresident(:,7)+Sresident(:,1), Sresident(:,8)+Sresident(:,2), Sresident(:,9)+Sresident(:,3), 'b', 'Marker', '+', 'MarkerIndices', 1:400:size(Sresident,1)); 
    t = plot3(Sgate(:,7)+Sgate(:,1), Sgate(:,8)+Sgate(:,2), Sgate(:,9)+Sgate(:,3), 'k', 'Marker', '*', 'MarkerIndices', 1:500:size(Sgate,1));
    g = plot3(Sc1(:,1)+Sc1(:,7), Sc1(:,2)+Sc1(:,8), Sc1(:,3)+Sc1(:,9), 'm');
    h = plot3(Sc2(:,1), Sc2(:,2), Sc2(:,3), 'm');
    plot3(St_transfer(:,1), St_transfer(:,2), St_transfer(:,3), 'm');
    plot3(St_departure(:,1), St_departure(:,2), St_departure(:,3), 'm');
    plot3(Stf(1:5000,1), Stf(1:5000,2), Stf(1:5000,3), 'm');
    plot3(St0(1:2000,1), St0(1:2000,2), St0(1:2000,3), 'm');
    scatter3(L(1,1), L(2,1), 0, 'k', 'filled');
    scatter3(L(1,2), L(2,2), 0, 'k', 'filled');
    scatter3(1-mu, 0, 0, 'k', 'filled');
    text(L(1,1)+1e-3, L(2,1)+1e-3, 5e-3, '$L_1$');
    text(L(1,2)+1e-3, L(2,2), 5e-3, '$L_2$');
    text(1-mu+1e-3, 0, 5e-3, '$M_2$');
    xlabel('Synodic $x$ coordinate');
    ylabel('Synodic $y$ coordinate');
    zlabel('Synodic $z$ coordinate');
    grid on;
    legend('Loiter orbit', 'Gateway orbit', 'Chaser trajectory', 'Location', 'northeast', 'AutoUpdate', 'off');
    title('Rendezvous simulation');

    for i = 1:dh:min([size(Sgate,1), size(Sresident,1)])
        T = scatter3(Sgate(i,1), Sgate(i,2), Sgate(i,3), 30, 'b', 'filled'); 
        V = scatter3(Sgate(i,1)+Sgate(i,7), Sgate(i,2)+Sgate(i,8), Sgate(i,3)+Sgate(i,9), 30, 'r', 'filled');
        H = scatter3(Sresident(i,1)+Sresident(i,7), Sresident(i,2)+Sresident(i,8), Sresident(i,3)+Sresident(i,9), 30, 'r', 'filled');
        J = scatter3(Sresident(i,1), Sresident(i,2), Sresident(i,3), 30, 'b', 'filled');

        drawnow;
        frame = getframe(W);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); 
        if (i == 1) 
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1e-3); 
        else 
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1e-3); 
        end 
        delete(T); 
        delete(V);
        delete(H)
        delete(J)
    end
    hold off
end