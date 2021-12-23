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

%Geostationary orbit 
theta = (0:1e-3:2*pi).';
Sgeo = hd*[cos(theta) sin(theta) zeros(length(theta),1)];

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
GNC.Control.MLQR.FloquetModes = diag(log(diag(lambda))/target_orbit.Period(end));
for i = 1:size(E,2)
    E(:,i) = E(:,i)/lambda(i,i);
end
GNC.Control.MLQR.FloquetDirections = E; 

%Initial conditions 
k = dimensionalizer(Lem, 1, 1, 10, 'Position', 0);      %Noise gain
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
Sresident = [Staux1; Staux2(2:end,:)];

%Control law
[~, ~, u] = GNCc_handler(GNC, Staux1(:,1:n), Staux1(:,n+1:end), tspan(tspan < 0.1));

%Control integrals
energy_sk = control_effort(tspan(tspan < 0.1), u, false);

total_cost(1) = 1000^1.025*norm(energy_sk(:,3));

%% Rendezvous with the resident
%Long-range rendezvous with MPC 
r_t0 = S(500,1:6);                            %Initial resident conditions
r_c0 = sHalo.';                               %Initial transfer vehicle conditions 
rho0 = r_c0-r_t0;                             %Initial relative conditions
s0 = [r_t0 rho0].';                           %Initial conditions of the target and the relative state

%Set up of the optimization
method = 'NPL';                               %Method to solve the problem
core = 'Linear';                              %Number of impulses
TOF = 0.1;                                    %Time of flight
cost_function = 'Position';                   %Target a position rendezvous

%Thruster characteristics 
Tmin = -1e-1;                                 %Minimum thrust capability (in velocity impulse)
Tmax = 1e-1;                                  %Maximum thrust capability (in velocity impulse)

%Main computation 
tspan = 0:dt:TOF; 
tic
[St_mpc, dV, state] = MPC_control(mu, cost_function, Tmin, Tmax, TOF, s0, core, method);
toc

%Control integrals
effort_mpc = control_effort(tspan, dV, true);

%Error in time 
[e_mpc, merit_mpc] = figures_merit(tspan, St_mpc);

%Regression of the resident stationkeeping orbit
tspan = 0:dt:pi/2-TOF; 
order = 50; 
[Cp, Cv, Cg] = CTR_guidance(order, tspan, Sresident(500:500+length(tspan)-1,7:12));

%Continuous close-rangee
GNC.Algorithms.Navigation = '';                             %Navigation algorithm

%Guidance parameters 
GNC.Algorithms.Control = 'SDRE';                            %Control algorithm
GNC.Algorithms.Guidance = 'CTR';               	            %Guidance algorithm
GNC.Guidance.CTR.Order = order;                             %Order of the approximation
GNC.Guidance.CTR.TOF = tspan(end);                          %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	        %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;                 %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;             %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.IntegralCoefficients = zeros(3,order);     %Coefficients of the Chebyshev approximation

GNC.System.mu = mu;                                         %System reduced gravitational parameter
GNC.System.Libration = [Ln gamma];                          %Libration point ID
GNC.Control.SDRE.Model = 'RLM';                             %SDRE model
GNC.Control.SDRE.Q = 2*eye(9);                              %Penalty on the state error
GNC.Control.SDRE.M = eye(3);                                %Penalty on the control effort

%Initial conditions 
s0 = [St_mpc(end,1:2*n) zeros(1,3)];                        %Initial conditions for the transfer vehicle

%Re-integrate trajectory
tic
[t, St] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Control law
[~, ~, u] = GNC_handler(GNC, St(:,1:6), St(:,7:15), t);   

%Control effort 
effort_sdre = control_effort(t, u, false);

%Total cost of the phase 
total_cost(2) = 1000*1.05*(norm(effort_sdre(:,3))+norm(effort_mpc(:,3)));

%Phase 1 trajectory 
S1c = [St_mpc(:,1:2*n); St(:,1:2*n)];       %Transfer vehicle trajectory
S1t = Sresident(500:size(S1c)-501,1:2*n);   %Resident spacecraft trajectory

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

%% Gateway stationkeeping orbit 
%Integration of the model
s0 = [target_orbit.Seeds(end,:) reshape(eye(n), [1 n^2])];
tspan = 0:dt:target_orbit.Period(end);
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);
Sn = repmat(S, floor(tf/target_orbit.Period(end))+1, 1);     

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
k = dimensionalizer(Lem, 1, 1, 10, 'Position', 0);      %Noise gain
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
[Shr, dVhr] = HTRC_guidance(mu, sequence, branch, rho, target_orbit, loiter_orbit, TOF, position_fixed, graphics);

%Control integrals
energy_hr = control_effort(NaN, dVhr, true);

%Guidance law
Sg = Shr.Trajectory; 

%Integration of the target trajectory
tspan = 0:dt:dt*(size(Sg,1)-1);
r_t0 = Sn(1,1:6);                                                                    %Initial target conditions
[~, St0] = ode113(@(t,s)cr3bp_equations(mu, 1, false, t, s), tspan, r_t0, options);  %Natural target trajectory

%Regression of the guidance coefficients 
Sgr = Sg-St0;                                                       %Relative guidance law
order = 300;                                                        %Order of the approximation
[Cp, Cv, Cg] = CTR_guidance(order, tspan, Sgr);

%Guidance parameters 
GNC.Guidance.CTR.Order = order;                             %Order of the approximation
GNC.Guidance.CTR.TOF = tspan(end);                          %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	        %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;                 %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;             %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.IntergralCoefficients = zeros(3,1);        %Coefficients of the Chebyshev approximation

%SMC tracking 
GNC.Algorithms.Guidance = '';               	            %Guidance algorithm
GNC.Algorithms.Navigation = '';                             %Navigation algorithm
GNC.Algorithms.Control = 'SMC';                             %Control algorithm
GNC.Guidance.Dimension = 9;                                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                                  %Dimension of the control law
GNC.System.mu = mu;                                         %System reduced gravitational parameters

GNC.Control.SMC.Parameters = [1.0000 0.6368 0.0008 0.0941];

%Re-integrate the trajectory
s0 = St_mpc(end,1:2*n);
tic 
[~, St_transfer] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Control law
[~, ~, u] = GNC_handler(GNC, St_transfer(:,1:6), St_transfer(:,7:12), tspan);    

%Control integrals
effort_smc = control_effort(tspan, u, false);

%Totalc cost of the maneuver 
total_cost(4) = 1000*1.025*(norm(effort_smc(:,3))+norm(energy_hr(:,3)));

%% Rendezvous with the Gateway 
%Compute a phasing arc using the CML guidance scheme
Ln = 2;                             %L2 ID
tf = 0.9*pi/2;                      %Interception time
Tsyn = target_orbit.Period(end);    %Synodic period
constraint.Flag = true;     
constraint.Period = Tsyn; 

initial_state = [Sgate(500,1:n)+Sgate(500,n+1:2*n) St_transfer(end,1:n)+St_transfer(end,n+1:2*n)];

[phasing_arc, dV, state_cml, S0] = CMC_guidance(mu, Ln, gamma, tf, constraint, initial_state, tol);

tspan = 0:dt:tf; 
effort_cml = control_effort(NaN, dV, true);

%Regression of the phasing arc
order = 50; 
[Cp, Cv, Cg] = CTR_guidance(order, tspan, phasing_arc);

%Guidance parameters 
GNC.Guidance.CTR.Order = order;                             %Order of the approximation
GNC.Guidance.CTR.TOF = tspan(end);                          %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	        %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;                 %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;             %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.IntergralCoefficients = zeros(3,order);    %Coefficients of the Chebyshev approximation

%SMC tracking 
GNC.Algorithms.Guidance = '';               	%Guidance algorithm
GNC.Algorithms.Navigation = '';                 %Navigation algorithm
GNC.Algorithms.Control = 'SMC';                 %Control algorithm
GNC.Guidance.Dimension = 9;                     %Dimension of the guidance law
GNC.Control.Dimension = 3;                      %Dimension of the control law
GNC.System.mu = mu;                             %System reduced gravitational parameters

GNC.Control.SMC.Parameters = [1.0000 0.6368 0.0008 0.0941];

%Re-integrate the trajectory
s0 = phasing_arc(end,1:2*n);
tic 
[~, S2_smc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
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
tf = 0.2;                                     %Time to dock                
times = [0 tf*rand(1,5)];                     %Times to impulse the spacecraft

%Compute the control law
impulses.Number = length(times);              %Number of impulses
impulses.Weights = eye(impulses.Number*3);    %Weightening matrix
impulses.Times = times;                       %Impulses times

cost = 'Position';                            %Cost function to target

%Controller scheme
tic
[S2_miss, dV, state] = MISS_control(mu, tf, S2_smc(end,1:2*n), tol, cost, impulses);
toc

%Control effort 
effort_miss = control_effort(NaN, dV, true);

%Totalc cost of the maneuver 
total_cost(5) = 1000*1.025*norm(effort_smc(:,3)+norm(effort_miss(:,3)));

%Complete rendezvous trajectory 
S2 = [S2_smc(:,1:2*n); S2_miss(:,1:2*n)];

%% Heteroclinic departure connection between the loitering and target orbit
%Manifold definition
TOF = pi;                       %Time of flight
rho = 90;                       %Manifold fibers to compute 

%Computation flags
position_fixed = true;          %Flag to determine a final target state
graphics = false;               %Flag to plot the manifolds

%Final target state
target_orbit.Trajectory = Sn;
target_orbit.TargetState = shiftdim(S2(end,1:n)); 

%Connection itenerary
sequence = [2 1];               %Connection itenerary
branch = ['R' 'R'];             %Manifold branches to be propagated

%Trajectory design core
[Shr, dVhr] = HTRC_guidance(mu, sequence, branch, rho, target_orbit, loiter_orbit, TOF, position_fixed, graphics);

%Control integrals
energy_hr = control_effort(NaN, dVhr, true);

%Guidance law
Sg = Shr.Trajectory; 

%Integration of the target trajectory
tspan = 0:dt:dt*(size(Sg,1)-1);
r_t0 = Sn(1,1:6);                                                                    %Initial target conditions
[~, St0] = ode113(@(t,s)cr3bp_equations(mu, 1, false, t, s), tspan, r_t0, options);  %Natural target trajectory

%Regression of the guidance coefficients 
Sgr = Sg-St0;                                                       %Relative guidance law
order = 300;                                                        %Order of the approximation
[Cp, Cv, Cg] = CTR_guidance(order, tspan, Sgr);

%Guidance parameters 
GNC.Guidance.CTR.Order = order;                             %Order of the approximation
GNC.Guidance.CTR.TOF = tspan(end);                          %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	        %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;                 %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;             %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.IntergralCoefficients = zeros(3,1);        %Coefficients of the Chebyshev approximation

%SMC tracking 
GNC.Algorithms.Guidance = '';               	            %Guidance algorithm
GNC.Algorithms.Navigation = '';                             %Navigation algorithm
GNC.Algorithms.Control = 'SMC';                             %Control algorithm
GNC.Guidance.Dimension = 9;                                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                                  %Dimension of the control law
GNC.System.mu = mu;                                         %System reduced gravitational parameters

GNC.Control.SMC.Parameters = [1.0000 0.6368 0.0008 0.0941];

%Re-integrate the trajectory
s0 = St_mpc(end,1:2*n);
tic 
[~, St_departure] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
toc

%Control law
[~, ~, u] = GNC_handler(GNC, St_departure(:,1:6), St_departure(:,7:12), tspan);    

%Control integrals
effort_smc = control_effort(tspan, u, false);

%Totalc cost of the maneuver 
total_cost(6) = 1000*1.025*(norm(effort_smc(:,3))+norm(energy_hr(:,3)));

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