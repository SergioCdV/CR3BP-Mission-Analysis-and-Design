%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 24/12/21 % 

%% Dynamical control Mission Scenario %% 
% This script provides an interface to test the RTMC to rendezvous two spacecraft and the Floquet mode strategy for collision avoidance maneuvers. 

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

%Time span 
dt = 1e-3;                          %Time step
tf = linspace(0.2, 0.8, 10);        %Rendezvous time

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
Az = 120e6;                                                 %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
T = target_orbit.Period;

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                    %Number of orbits to continuate
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, T};                           %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
direction = 1;                                              %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                         %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(2,:), maxIter, tol);

%% Modelling in the synodic frame %%
r_t0 = target_orbit.Trajectory(500,1:6);        %Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);          %Initial chaser conditions 
rho0 = r_c0-r_t0;                               %Initial relative conditions
s0 = [r_t0 rho0].';                             %Initial conditions of the target and the relative state

%Integration of the model
tspan = 0:dt:2*pi;
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                      %Reconstructed chaser motion via Encke method

%% GNC: TITA %%
%Differential corrector set up
tol = 1e-10;                                    %Differential corrector tolerance

%Cost function matrices
penalties.R = eye(3);                           %Penalty on the impulse
penalties.Q = 2*eye(6);                         %Penalty on the state error
penalties.M  = 0.1*eye(6);                      %Penalty on the state noise

%Select measuring times 
target_points.Noise = true;                     %Boolean to account for state noise
thruster_model.Sigma = 0.01;                    %Velocity noise dependance on the velocity impulse

%Rotational misalignment of the thrusters
thruster_model.Rotation = [1 0 0; 0 cos(pi/18) sin(pi/18); 0 -sin(pi/18) cos(pi/18)];           

%Cost function 
cost_function = 'Position';                     %Cost function to target
two_impulsive = true;                           %Two-impulsive rendezvous boolean

%Preallocation 
V = zeros(1,length(tf));                        %Total maneouvre cost

for i = 1:length(tf)
    tspan = 0:dt:tf(i);                         %Integration time span
    target_points.Times = tf(i)*rand(1,3);      %Times to measure the state noise

    tic
    [St, dV, state] = TITA_control(mu, tf(i), s0, cost_function, zeros(1,3), two_impulsive, ...
                                   penalties, target_points, thruster_model, tol);
    toc

    %Error performance 
    [e, merit] = figures_merit(tspan, St);

    %Total maneuver metrics 
    effort_rtps = control_effort(NaN, dV, true);

    V(i) = norm(effort_rtps(:,3))*1.05*1000;
end
                           
%% Results %% 
%Plot results 
figure(1) 
view(3) 
hold on 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'b'); 
plot3(S(:,1), S(:,2), S(:,3), 'r'); 
r = plot3(St(:,1)+St(:,7), St(:,2)+St(:,8), St(:,3)+St(:,9), 'm');
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
text(L(1,Ln)+1e-3, L(2,Ln), 5e-3, '$L_2$');
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Initial orbit', 'Target orbit', 'Rendezvous arc', 'Location', 'northwest');

%Rendezvous animation 
if (false)
    figure(5) 
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