%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 3: LQR/SDRE control law %% 
% This script provides an interface to test the two-impusilve rendezvous strategia using the STM of 
% the dynamical model. 

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% The classical LQR for a time-varying model is solved to drive the relative phase space vector to the origin 
% (rendezvous condition).

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
mass = 1e-21;

%Time span 
dt = 1e-3;                          %Time step
tf = 0.6;                           %Rendezvous time
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

%% Modelling in the synodic frame %%
index = fix(tf/dt);                                         %Rendezvous point
if (index > size(target_orbit.Trajectory,1))
    index = mod(index, size(target_orbit.Trajectory,1));    %Rendezvous point
end
r_t0 = target_orbit.Trajectory(33,1:6);                     %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% GNC: LQR/SDRE control law %%
%Approximation 
order = 2;                                                  %Order of the approximation 

%Cost function matrices 
Q = eye(n);                                                 %Cost weight on the state error
R = eye(n/2);                                               %Cost weight on the spent energy

%Preallocation 
S = zeros(length(tspan), 2*n);                              %Relative orbit trajectory

%Input matrix 
B = (1/mass)*[zeros(n/2); eye(n/2)];

%Compute the trajectory
for i = 1:length(tspan)
    %Compute the relative Legendre coefficient c2 
    rc = relegendre_coefficients(mu, S(i,1:3).', order); 
    c2 = rc(2); 
    
    %Evaluate the linear model 
    A = [zeros(3) eye(3); 1+2*c2 0 0 0 2 0; 0 1-c2 0 -2 0 0; 0 0 0 -c2 0 0];     %Linear model state matrix
    
    %Compute the LQR matrix 
    [~, ~, K] = care(A, B, Q, R);                                                %Feedback gain matrix
    
    %Compute the control law
    u = -K*shiftdim(S(i,7:12));                                                  %Feedback control law
    
    %Re-integrate trajectory
    [~, s] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke C', t, s, u), [0 dt], S(i,:), options);
    
    %Update initial conditions
    S(i+1,:) = s(end,:);
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
plot3(S(:,7), S(:,8), S(:,9)); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative motion in the configuration space');

%Rendezvous animation 
if (false)
    figure(3) 
    view(3) 
    grid on;
    hold on
    plot3(Sn(1:index,1), Sn(1:index,2), Sn(1:index,3), 'k-.'); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    title('Rendezvous simulation');
    for i = 1:size(St,1)
        T = scatter3(Sn(i,1), Sn(i,2), Sn(i,3), 30, 'b'); 
        V = scatter3(St(i,1)+St(i,7), St(i,2)+St(i,8), St(i,3)+St(i,9), 30, 'r');

        drawnow;
        delete(T); 
        delete(V);
    end
hold off
end

