%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 25/08/22 % 

%% Center Manifold Lissajous Guidance demonstration for IAC 2022 %% 
% This script provides an interface to demonstrate the CML guidance core.

% The relative motion of two spacecraft in the same halo orbit around L1 in the
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
tf = 1;                             %Rendezvous time
tspan = 0:dt:pi;                    %Integration time span

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
Az = 20e6;                                                          %Orbit amplitude out of the synodic plane. 
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
num = 5;                                                            %Number of orbits to continuate
method = 'SPC';                                                     %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                        %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};                 %Object and characteristics to continuate
corrector = 'Plane Symmetric';                                      %Differential corrector method
direction = 1;                                                      %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                                 %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);

%Halo characteristics 
Az = 10e6;                                                          %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Modelling in the synodic frame %%
r_t0 = target_orbit.Trajectory(1,1:6);                            %Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);                              %Initial chaser conditions 
rho0 = r_c0-r_t0;                                                   %Initial relative conditions
s0 = [r_t0 rho0].';                                                 %Initial conditions of the target and the relative state

%Integration of the model
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
Sn = S;                

%Reconstructed chaser motion 
Sr = S(:,1:6)+S(:,7:12);                                            %Reconstructed chaser motion via Encke method

%% Generate the guidance trajectory
%CML Guidance trajectory
Tsyn = target_orbit.Period*chaser_orbit.Period/(target_orbit.Period+chaser_orbit.Period);
constraint.Flag = false; 
constraint.Period = Tsyn; 

[Str, V1, ilg_state, S0] = CMLG_guidance(mu, Ln, gamma, tf, constraint, [r_t0 r_c0], tol);
Str(end,10:12) = zeros(1,3);

% Error evaluation 
[error, merit] = figures_merit(tspan, Str);

% Comparison against the TISS controller
[Str2, dV, tiss_state] = TISS_control(mu, tf, [r_t0 r_c0-r_t0].', tol, 'Position', true);  

% Final absolute trajectories
tspan = 0:dt:tf;                  % Final transfer time span
S0 = S0(:,1:6)+S0(:,7:12);          % Natural quasi-periodic model    
St = Str(:,1:n)+Str(:,n+1:2*n);     % Transfer trajectory

% Error evaluation 
[error, merit2] = figures_merit(tspan, Str2);

%% Results 
%Plot results 
figure(1) 
view(3) 
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'b', 'LineWidth', 0.9); 
plot3(Sr(:,1), Sr(:,2), Sr(:,3), '-ob', 'LineWidth', 0.9, 'MarkerIndices', floor(linspace(1,size(Sr,1),5))); 
plot3(St(:,1), St(:,2), St(:,3), 'r', 'LineWidth', 1); 
plot3(S0(:,1), S0(:,2), S0(:,3), 'g', 'LineWidth', 0.5); 
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;
legend('Reference target orbit', 'Chaser orbit', 'Guidance transfer orbit', 'Quasi-periodic guess', 'AutoUpdate', 'off')
plot3(L(1,Ln), L(2,Ln), 0, '+k');
labels = {'$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$'};
text(L(1,Ln)-1e-3, L(2,Ln)-1e-3, 1e-2, labels{Ln});
hold off

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan(1:size(Str,1)), Str(:,7)); 
plot(tspan(1:size(Str,1)), Str(:,8)); 
plot(tspan(1:size(Str,1)), Str(:,9)); 
hold off
xlabel('$t$');
ylabel('$\mathbf{\rho}$');
grid on;
legend('$x$', '$y$', '$z$');
subplot(1,2,2)
hold on
plot(tspan(1:size(Str,1)), Str(:,10)); 
plot(tspan(1:size(Str,1)), Str(:,11)); 
plot(tspan(1:size(Str,1)), Str(:,12)); 
hold off
xlabel('$t$');
ylabel('$\dot{\mathbf{\rho}}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');

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
        drawnow;
        delete(T); 
        delete(V);
    end
    hold off
end
