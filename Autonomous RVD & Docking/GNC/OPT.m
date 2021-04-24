%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 7: MPC guidance-control law %% 
% This script provides an interface to test MPC control law rendezvous strategies for
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

%Spacecraft mass 
mass = 1e-10;

%Time span 
dt = 1;                             %Time step
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
index = fix(tf/dt);                                         %Rendezvous point
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% Optimal control problem -equivalent to MPC- 
%Simulation time bounds
P.bounds.initialTime.low = 0;
P.bounds.initialTime.upp = 0;
P.bounds.finalTime.low = 0;
P.bounds.finalTime.upp = tf;

%Boundary conditions bounds
P.bounds.state.low = -inf*ones(6,1);
P.bounds.state.upp = inf*ones(6,1);
P.bounds.initialState.low = rho0.';
P.bounds.initialState.upp = rho0.';
P.bounds.finalState.low = zeros(6,1);
P.bounds.finalState.upp = zeros(6,1);

%Control law  bounds
P.bounds.control.low = [0; 0; 0];
P.bounds.control.upp = 0.1*ones(3,1);

%Initial guess
P.guess.time = [0 tf];  
P.guess.state = [rho0.' zeros(6,1)];
P.guess.control = [P.bounds.control.low P.bounds.control.low];

%Dynamics function
P.func.dynamics = @(t,x,u)(opt_model(mu, Sn, t, x, u));   

%Boundary constraint
P.func.bndCst = @(t0,x0,tF,xF)(pathConstraint(xF));

%Objective function
tol = 1e-3;
P.func.pathObj = @(t,x,u)(dot(u,u,1));

%Select transcription method
method = 'rungeKutta';
switch (method)
    case 'trapezoid'
        P.options(1).method = 'trapezoid';
        P.options(1).defaultAccuracy = 'medium';
        P.options(1).nlpOpt.MaxFunEvals = 2e5;
        P.options(1).nlpOpt.MaxIter = 1e5;
        P.options(2).method = 'trapezoid';
        P.options(2).defaultAccuracy = 'medium';
        P.options(2).nlpOpt.MaxFunEvals = 2e4;
        P.options(2).nlpOpt.MaxIter = 1e5;
        
    case 'rungeKutta'
        P.options(1).method = 'rungeKutta';
        P.options(1).defaultAccuracy = 'low';
        P.options(1).rungeKutta.nSegment = 20;
        P.options(2).method = 'rungeKutta';
        P.options(2).defaultAccuracy = 'high';
        
    case 'chebyshev'
        P.options(1).method = 'chebyshev';
        P.options(1).defaultAccuracy = 'low';
        P.options(2).method = 'chebyshev';
        P.options(2).defaultAccuracy = 'low';
        P.options(2).chebyshev.nColPts = 15;
    otherwise 
        error('No valid transcription method was selected');
end

%Solve the problem
soln = optimTraj(P);

%Solution interpolation
t = linspace(soln(end).grid.time(1),soln(end).grid.time(end),250);
x = soln(end).interp.state(t);
u = soln(end).interp.control(t);

%% Auxiliary functions 
function [c, ceq] = pathConstraint(xF)
    c = []; 
    ceq = xF;
end