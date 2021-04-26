%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 20/02/21 % 

%% Nonlinear models %% 
% This script provides a script to test different nonlinear models for relative motion in the CR3BP. 

% The relative motion of two spacecraft in a halo orbit around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion % 

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
%Time span 
dt = 1e-3;                          %Time step
tmax = 0.8*pi;                      %Maximum time of integration (corresponding to a synodic period)
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
[halo_orbit, state] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Modelling in the synodic frame %% 
r_t0 = halo_orbit.Trajectory(1,1:6);                        %Initial target conditions
r_c0 = halo_orbit.Trajectory(2,1:6);                        %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0];                                           %Initial conditions of the target and the relative state

%Integration of the model
[~, S_c] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, r_c0, options);
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke', t, s), tspan, s0, options);
[t, Sn] = ode113(@(t,s)nlr_model(mu, true, false, 'Second order', t, s), tspan, s0, options);

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method
error = S_c-S_rc;                                           %Error via the full nonlinear model with the Encke method
S_rcn = Sn(:,1:6)+Sn(:,7:12);                               %Reconstructed chaser motion via the full nonlinear model
error_n = S_c-S_rcn;                                        %Error via the full nonlinear model

ep = zeros(size(error,1), 1);                               %Position error (L2 norm) via Encke's method
epn = zeros(size(error,1), 1);                              %Position error (L2 norm) via direct integration
for i = 1:size(error,1)
    ep(i) = norm(error(i,:));
    epn(i) = norm(error_n(i,:));
end

%% Results in the inertial frame %% 
%Preallocation 
inertial = zeros(size(S_c,1),3);      %Position error in the inertial frame

%Main computation
for i = 1:size(S_c,1)
     elaps = dt*(i-1);                                                 %Elaps time
     inertial(i,:) = (inertial2synodic(elaps, Sn(i,7:9).', 1)).';      %Relative position in the synodic frame
end

%% Results in the libration synodic frame %% 
%Preallocation 
libration = zeros(size(S_c,1),3);                                               %Position error in the synodic frame

%Main computation
for i = 1:size(S_c, 1)
     libration(i,:) = (synodic2lagrange(mu, gamma, Ln, Sn(i,7:9).', 1)).';      %Relative position in the synodic frame
end

%% Results in the Frenet-Serret frame of the target orbit %% 
%Preallocation
T = zeros(size(S_c,1),3,3);             %Synodic to Frenet frame rotation matrix
frenet = zeros(size(S_c,1), 3);         %Position error in the Frenet frame

%Main computation
for i = 1:size(S_c,1)
    T(i,1:3,1:3) = frenet_triad(mu, S(i,1:6));                %Frenet-Serret frame
    frenet(i,:) = (shiftdim(T(i,:,:)).'*S(i,7:9).').';        %Relative position in the Frenet-Serret frame of the target orbit
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
    U = -(mup(1)*((1/norm(rho+r_t-R(:,1)))-(dot(rho,(R(:,1)-r_t))/norm(r_t-R(:,1)))) ...
         +mup(2)*((1/norm(rho+r_t-R(:,2)))-(dot(rho,(R(:,2)-r_t))/norm(r_t-R(:,2)))));
     
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
plot3(S_rcn(:,1), S_rcn(:,2), S_rcn(:,3)); 
hold off
legend('Target motion', 'Chaser motion', 'New chaser motion'); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Reconstruction of the chaser motion');

figure(2) 
hold on
plot(t, log(ep), 'b'); 
plot(t, log(epn), 'r');
hold off
grid on
xlabel('Nondimensional epoch'); 
ylabel('Relative error in the synodic frame');
title('Error in the phase space vector (L2 norm)')
legend('Encke method error', 'Direct integration error');

figure(3) 
plot(t, H, 'b'); 
grid on
xlabel('Nondimensional epoch'); 
ylabel('Relative Hamiltonian');
title('Relative energy evolution')

%Relative orbit plots
if (false)
    figure(4) 
    plot3(inertial(:,1), inertial(:,2), inertial(:,3)); 
    xlabel('Nondimensional x coordinate'); 
    ylabel('Nondimensional y coordinate');
    zlabel('Nondimensional z coordinate');
    grid on
    title('Relative orbit in the inertial frame'); 

    figure(5) 
    plot3(libration(:,1), libration(:,2), libration(:,3)); 
    xlabel('Nondimensional x coordinate'); 
    ylabel('Nondimensional y coordinate');
    zlabel('Nondimensional z coordinate');
    grid on
    title('Relative orbit in the libration synodic frame'); 

    figure(6) 
    plot3(frenet(:,1), frenet(:,2), frenet(:,3)); 
    xlabel('Nondimensional x coordinate'); 
    ylabel('Nondimensional y coordinate');
    zlabel('Nondimensional z coordinate');
    grid on
    title('Relative orbit in the Frenet-Serret frame');
end
