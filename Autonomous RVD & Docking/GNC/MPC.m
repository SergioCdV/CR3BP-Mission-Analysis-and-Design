%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 6: MPC guidance-control law %% 
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

%% MPC guidance scheme
%Set up of the optimization
method = 'NPL';
model = 'Fixed libration';

%Environment characteristics 
cn = legendre_coefficients(mu, Ln, gamma, 2);   %Legendre coefficients

%Thruster characteristics 
Tmax = 0.1;

%Preallocation 
Sc = zeros(length(tspan),length(s0));           %Preallocate the trajectory
u = zeros(3,length(tspan));                     %Preallocate the control vector
e = zeros(1,length(tspan));                     %Preallocate the error vector

Sc(1,:) = s0;                                   %Initial conditions

for i = 1:1
    %Shrink the horizon 
    Dt = tspan(i:end); 
    
    %Compute the control law 
    U = MPC_guidance(method, model, mu, Sc(i,:), Dt, cn, Tmax);
%     u(:,i) = U(:,1); 
% 
%     %Re-integrate the trajectory 
%     [~, s] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke C', t, s, u(:,i)), [0 dt], Sc(i,:), options);
% 
%     %Update initial conditions
%     Sc(i+1,:) = s(end,:);
%     
%     %Update the error 
%     e(i) = norm(s(end,:));
end

%Shrink the integrated trajectory
Sc = Sc(1:end-1,:);

%% Results %% 
%Plot results 
figure(1) 
view(3) 
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3)); 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
hold off
legend('Target motion', 'Chaser motion'); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Reconstruction of the natural chaser motion');

%Plot relative phase trajectory
figure(2) 
view(3) 
plot3(Sc(:,7), Sc(:,8), Sc(:,9)); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative motion in the configuration space');

%Configuration space error 
figure(3)
plot(tspan, e); 
xlabel('Nondimensional epoch');
ylabel('Absolute error');
grid on;
title('Absolute error in the configuration space (L2 norm)');

%Rendezvous animation 
if (false)
    figure(4) 
    view(3) 
    grid on;
    hold on
    plot3(Sc(1:index,1), Sc(1:index,2), Sc(1:index,3), 'k-.'); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    title('Rendezvous simulation');
    for i = 1:size(Sc,1)
        T = scatter3(Sc(i,1), Sc(i,2), Sc(i,3), 30, 'b'); 
        V = scatter3(Sc(i,1)+Sc(i,7), Sc(i,2)+Sc(i,8), Sc(i,3)+Sc(i,9), 30, 'r');

        drawnow;
        delete(T); 
        delete(V);
    end
hold off
end

%% Auxiliary functions
%Controlability analysis function 
function [commands] = MPC_guidance(method, model, mu, s0, tspan, cn, Tmax)
    %Constants 
    finalHorizonIndex = length(tspan);
    
    %Cost function 
    cosfunc = @(u)(-norm(u(1:3,:)));
    
    %Linear constraints 
    A = []; 
    b = []; 
    Aeq = []; 
    beq = [];
    
    switch (method)
        case 'Genetic algorithm'
            dof = length(tspan);
            PopSize = 100;          %Population size for each generation
            MaxGenerations = 10;    %Maximum number of generations for the evolutionary algorithm
            options = optimoptions(@ga,'PopulationSize', PopSize, 'MaxGenerations', MaxGenerations, 'ConstraintTolerance', 1e-1, 'PlotFcn', @gaplotbestf);
            lb = zeros(3*finalHorizonIndex,1);
            ub = Tmax*ones(3*finalHorizonIndex,1);
            commands = ga(@(u)cosfunc(u), dof, A, b, Aeq, beq, lb, ub, @(u)nonlcon(model, mu, s0, tspan, cn, u), options);
        case 'NPL'
            lb = zeros(3,finalHorizonIndex);
            ub = Tmax*ones(3,finalHorizonIndex);
            u0 = Tmax*ones(3,length(tspan));
            commands = fmincon(@(u)cosfunc(u), u0, A, b, Aeq, beq, lb, ub, @(u)nonlcon(model, mu, s0, tspan, cn, u));
        otherwise 
            error('No valid method was chosen');
    end
end

function [c, ceq] = nonlcon(model, mu, s0, tspan, cn, u)
    %Approximation 
    n = 6;                              %Dimension of the state vector
    
    %Set up 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
    
    %Preallocation 
    S = zeros(length(tspan), 2*n);      %Preallocation of the trajectory
    S(1,:) = s0;
    
    %Integration of the trajectory 
    dt = tspan(2)-tspan(1); 
    for i = 1:length(tspan)
        [~,s] = ode113(@(t,s)lr_model(mu, cn, 1, false, model, t, s, true, u(:,i)), [0 dt], S(i,:), options);
        S(i+1,:) = s(end,:);
    end
    
    %Nonlinear constraints
    c = [];                         %Empty nonlinear inequalities
    ceq = S(end,n+1:2*n);           %Rendezvous conditions
end
