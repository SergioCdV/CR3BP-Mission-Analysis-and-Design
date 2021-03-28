%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 23/03/21 % 

%% Phase space study %% 
% This script provides an interface to test how the phase space volumens evolve. 

% The relative motion of two spacecraft in two halo orbits around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion 

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

%Continuate the firts halo to locate the chaser spacecraft
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                    %Number of orbits to continuate
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};         %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
direction = 1;                                              %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                         %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(end,1:6), maxIter, tol);

%% Modelling in the synodic frame %% 
r_t0 = target_orbit.Trajectory(1,1:6).';                    %Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6).';                    %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0; rho0];                                          %Initial conditions of the target and the relative state
ds0 = eye(length(r_t0));                                    %Initial phase space error vector
ds0 = reshape(ds0, [length(r_t0)^2 1]);                     %Initial phase space error vector
s0 = [s0; ds0];                                             %Initial conditions of the target and the relative state and the error

%Integration of the model
[~, S_c] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, r_c0, options);
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke V', t, s), tspan, s0, options);

%% Analysis of the state transition matrix evolution
%Preallocation
STM = zeros(size(S,1), length(r_t0), length(r_t0));         %STM of the system
detSTM = zeros(size(S,1),1);                                %Determinant of the STM
drho = zeros(size(S,1),3);                                  %Total variation along each configuration space coordinate
dv = zeros(size(S,1),3);                                    %Total variation along each velocity space coordinate

%Main computation
for i = 1:size(S,1)
    STM(i,:,:) = reshape(S(i,13:end), [length(r_t0), length(r_t0)]);
    detSTM(i) = det(shiftdim(STM(i,:,:)));
    drho(i,1:3) = [sum(shiftdim(STM(i,1,:))) sum(shiftdim(STM(i,2,:))) sum(shiftdim(STM(i,3,:)))];
    dv(i,1:3) = [sum(shiftdim(STM(i,4,:))) sum(shiftdim(STM(i,5,:))) sum(shiftdim(STM(i,6,:)))];
end

%% Analysis of the phase space volume invariancy 
%Preallocation
J = zeros(size(S,1), length(r_t0), length(r_t0));           %Jacobian all along the trajectory
detJ = zeros(size(S,1),1);                                  %Determinant of the Jacobian

%Main computation
for i = 1:size(S,1)
    J(i,:,:) = rel_jacobian(mu, S(i,:).');                  %Jacobian of the system
    detJ(i) = det(shiftdim(J(i,:,:)));                      %Determinant of the Jacobian
end

%Fourier transform of the determinant 
Adet = fft(detJ(1:11000));
f = (1/(dt*length(Adet)))*(1:length(Adet));
PSD = Adet.*conj(Adet)/length(Adet);

%% Results %% 
% Plot results 
figure(1) 
plot(t, log(drho(:,1:3))); 
title('Evolution of the position initial displacements'); 
grid on; 
legend('x error', 'y error', 'z error'); 
xlabel('Nondimensional time'); 
ylabel('Error evolution');

figure(2) 
plot(t, log(dv(:,1:3))); 
title('Evolution of the velocity initial displacements'); 
grid on; 
legend('x error', 'y error', 'z error');
xlabel('Nondimensional time'); 
ylabel('Error evolution');

figure(3)
plot(t, detSTM); 
title('Evolution of the determinant of the STM'); 
grid on; 
xlabel('Nondimensional time'); 
ylabel('Determinant of the STM');

figure(4)
plot(t, detJ); 
title('Evolution of the determinant of the Jacobian'); 
grid on; 
xlabel('Nondimensional time'); 
ylabel('Determinant of the Jacobian');

figure(5)
plot(f, PSD); 
title('Power spectral density of the determinant of the Jacobian'); 
grid on; 
xlabel('Nondimensional frequency'); 
ylabel('PSD');

if (false)
    figure(1) 
    view(3) 
    plot3(S(:,7), S(:,8), S(:,9), 'g');
    grid on;
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    title('Phase space evolution');
    hold on
    for i = 1:size(S,1) 
       X = quiver3(S(i,7), S(i,8), S(i,9), 1e-3*STM(i,1,1), 1e-3*STM(i,2,1), 1e-3*STM(i,3,1), 'r');
       Y = quiver3(S(i,7), S(i,8), S(i,9), 1e-3*STM(i,1,2), 1e-3*STM(i,2,2), 1e-3*STM(i,3,2), 'b');
       Z = quiver3(S(i,7), S(i,8), S(i,9), 1e-3*STM(i,1,3), 1e-3*STM(i,2,3), 1e-3*STM(i,3,3), 'k');
       drawnow;
       delete(X); 
       delete(Y);
       delete(Z);
    end
    hold off
end