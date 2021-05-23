%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 17/04/21 % 

%% GNC 5: SMC control law %% 
% This script provides an interface to test the Sliding Mode Control strategy for rendezvous missions. 

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% The SMC controller is solved to drive the relative phase space vector to the origin (rendezvous condition).

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
optimization = true;                %Optimize the controller parameters 

%Time span 
dt = 1e-3;                          %Time step
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
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(2,:), maxIter, tol);

%% Modelling in the synodic frame %%
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                  %Reconstructed chaser motion via Encke method

%% GNC algorithms definition 
GNC.Algorithms.Guidance = '';               %Guidance algorithm
GNC.Algorithms.Navigation = '';             %Navigation algorithm
GNC.Algorithms.Control = 'SMC';             %Control algorithm

GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
GNC.Control.Dimension = 3;                  %Dimension of the control law

GNC.System.mu = mu;                         %System reduced gravitational parameter

GNC.Control.SMC.Parameters = [1 1 0.9 0.1]; %Controller parameters

%% GNC: SMC control law %%
% %Preallocation 
e = zeros(1,length(tspan));             %Error to rendezvous 
energy = zeros(3,2);                    %Energy vector preallocation

%Re-integrate trajectory
[~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);

%Error in time 
for i = 1:length(tspan)
    e(i) = norm(Sc(i,7:12));
end

%Compute the error figures of merit 
ISE = trapz(tspan, e.^2);                                  %Integral of the square error                             
IAE = trapz(tspan, abs(e));                                %Integral of the absolute value of the error

%Compute the control effort
GNC.Control.SMC.Parameters = GNC.Control.SMC.Parameters;   %Controller parameters

%Control law
[~, ~, u] = GNC_handler(GNC, Sc(:,1:6), Sc(:,7:12));    
   
%Control integrals
for i = 1:size(u,1)
    energy(i,1) = trapz(tspan, u(i,:).^2);                 %L2 integral of the control
    energy(i,2) = trapz(tspan, sum(abs(u(i,:)),1));        %L1 integral of the control
end

%% Optimize the controller parameters
if (optimization)
    GNC.Control.SMC.Parameters = [1 SMC_optimization(mu, 'L2', s0, tf)];

    %Preallocation 
    e2 = zeros(1,length(tspan));             %Error to rendezvous 
    energy2 = zeros(3,2);                    %Energy vector preallocation

    %Re-integrate the trajectory
    [~, Sc2] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);

    %Error in time 
    for i = 1:length(tspan)
        e2(i) = norm(Sc2(i,7:12));
    end

    %Compute the error figures of merit 
    ISE2 = trapz(tspan, e2.^2);                                  %Integral of the square error                             
    IAE2 = trapz(tspan, abs(e2));                                %Integral of the absolute value of the error

    %Control law
    [~, ~, u2] = GNC_handler(GNC, Sc2(:,1:6), Sc2(:,7:12));    

    %Control integrals
    for i = 1:size(u,1)
        energy2(i,1) = trapz(tspan, u2(i,:).^2);                 %L2 integral of the control
        energy2(i,2) = trapz(tspan, sum(abs(u2(i,:)),1));        %L1 integral of the control
    end

    %Save results for statistics purposes 
    fileID = fopen('Autonomous RVD & Docking\Mission scenarios\smc_parameters.txt', 'a+'); 
    output = [GNC.Control.SMC.Parameters Sc(end,7:12) Sc2(end,7:12) ISE IAE ISE2 IAE2 ...
              norm(energy(:,1)) norm(energy(:,2)) norm(energy2(:,1)) norm(energy2(:,2))];
    format = '%.6f %.6f %.6f %.6f && %.6f %.6f %.6f %.6f %.6f %.6f && %.6f %.6f %.6f %.6f %.6f %.6f && %.6f %.6f && %.6f %.6f && %.6f %.6f && %.6f %.6f \n';
    fprintf(fileID, format, output);
    fclose(fileID); 
end

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

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan, Sc(:,7)); 
plot(tspan, Sc(:,8)); 
plot(tspan, Sc(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('x coordinate', 'y coordinate', 'z coordinate');
title('Relative position evolution');
subplot(1,2,2)
hold on
plot(tspan, Sc(:,10)); 
plot(tspan, Sc(:,11)); 
plot(tspan, Sc(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinate');
grid on;
legend('x velocity', 'y velocity', 'z velocity');
title('Relative velocity evolution');

%Configuration space evolution
figure(6)
subplot(1,2,1)
hold on
plot(tspan, Sc2(:,7)); 
plot(tspan, Sc2(:,8)); 
plot(tspan, Sc2(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('x coordinate', 'y coordinate', 'z coordinate');
title('Optimal relative position evolution');
subplot(1,2,2)
hold on
plot(tspan, Sc2(:,10)); 
plot(tspan, Sc2(:,11)); 
plot(tspan, Sc2(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Optimal relative velocity coordinate');
grid on;
legend('x velocity', 'y velocity', 'z velocity');
title('Relative velocity evolution');

%Configuration space error 
figure(4)
hold on
plot(tspan, log(e)); 
plot(tspan, log(e2)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Absolute error');
grid on;
title('Absolute error in the configuration space (L2 norm)');
legend('Nominal error', 'Optimal error')

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