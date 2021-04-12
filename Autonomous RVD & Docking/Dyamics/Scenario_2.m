%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 24/03/21 % 

%% Scenario 2 %% 
% This script provides an interface to test different nonlinear models for relative motion in the CR3BP
% in different orbits. 

% The relative motion of two spacecraft in two halo orbits around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion % 

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

% Relative motion happens between a halo orbit for the target and random
% initial conditions (pseudo-halo) for the chaser. 

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
%Time span 
dt = 1e-3;                          %Time step
tmax = 4*pi;                        %Maximum time of integration (corresponding to a synodic period)
tspan = 0:dt:tmax;                  %Integration time span

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
Az = 200e6;                                                 %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                    %Number of orbits to continuate
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};         %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
direction = 1;                                              %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                         %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);

%% Modelling in the synodic frame %% 
r_t0 = target_orbit.Trajectory(1,1:6);                      %Initial target conditions
r_c0 = chaser_seed.Seeds(1,1:6);                            %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0];                                           %Initial conditions of the target and the relative state

%Integration of the model
[~, S_c] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, r_c0, options);
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke', t, s), tspan, s0, options);

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%Error in the computation
error = S_c-S_rc;                                           %Position error (L2 norm) via Encke's method          
ep = zeros(size(error,1), 1);                               
for i = 1:size(error,1)
    ep(i) = norm(error(i,:));
end

%% Evolution of the Hamiltonian of the system
%Preallocation
H = zeros(size(S,1),1);               %Relative Hamiltonian

%Main computation
for i = 1:size(S,1)
    r_t = S(i,1:3).';                 %Target position       
    rho = S(i,7:9).';                 %Relative position
    v = S(i,10:12).';                 %Relative velocity
    R = [-mu 1-mu; 0 0; 0 0];         %Position of the primaries
    mup = [1-mu; mu];                 %Reduced gravitational parameters of the primaries
    
    %Kinetic energy
    T = (1/2)*norm(v-[0 -1 0; 1 0 0; 0 0 0]*rho)^2;
    
    %Potential energy
    U = 0;
    for j = 1:length(mup)
        U = U-mup(j)*((1/norm(rho+r_t-R(:,j)))-(dot(rho,(R(:,j)-r_t))/norm(r_t-R(:,j))^3));
    end
     
    %Relative Hamiltonian
    H(i) = T+U;
end

%% Results %% 
% Plot results 
figure(1) 
view(3) 
hold on
plot3(S(:,1), S(:,2), S(:,3)); 
plot3(S_c(:,1), S_c(:,2), S_c(:,3)); 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
hold off
legend('Target motion', 'Chaser motion', 'New chaser motion'); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Reconstruction of the chaser motion');

figure(2) 
view(3) 
hold on
plot3(S(:,7), S(:,8), S(:,9));  
hold off 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative orbit');

figure(3) 
hold on
plot(t, log(ep), 'b'); 
hold off
grid on
xlabel('Nondimensional epoch'); 
ylabel('Relative error in the synodic frame');
title('Error in the phase space vector (L2 norm)')
legend('Encke method error');

figure(4) 
plot(t, H, 'b'); 
grid on
xlabel('Nondimensional epoch'); 
ylabel('Relative Hamiltonian');
title('Relative energy evolution')