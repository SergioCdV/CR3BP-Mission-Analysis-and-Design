%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 1: Two-impulsive rendezvous via STM %% 
% This script provides an interface to test the two-impusilve rendezvous strategy using the STM of 
% the dynamical model. 

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% The first impulse is known as targetting, aimed to nullify the relative
% position between chaser and target after some flight time tf. The second
% impulse nullifies the relative velocity between the two to complete the
% rendezvous.

% In the relative phase space, the relative particle is driven to the
% origin of the synodic frame.

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
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state
Phi = eye(length(r_t0));                                    %Initial STM
Phi = reshape(Phi, [length(r_t0)^2 1]);                     %Initial STM
s0 = [s0; Phi];                                             %Initial conditions

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke V', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% GNC: two impulsive rendezvous, single shooting scheme %%
%Differential corrector set up
S = S(1:index,:);                   %Restrict the time integration span
maxIter = 20;                       %Maximum number of iterations
tol = 1e-10;                        %Differential corrector tolerance
GoOn = true;                        %Convergence boolean 
iter = 1;                           %Initial iteration 

%Preallocation 
dV = zeros(3,maxIter);              %Targeting impulse

%First impulse: targeting 
while ((GoOn) && (iter < maxIter))
    %Initial impulse
    error = S(end,7:9).';                                           %Velocity error
    STM = reshape(S(end,13:end), [length(r_t0) length(r_t0)]);      %STM evaluated at time tf
    dV(:,iter) = STM(1:3,4:6)\error;                                %Required impulse
    
    %Dynamics integration 
    s0(10:12) = s0(10:12)-dV(:,iter);                                                       %New initial conditions
    [~, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke V', t, s), tspan, s0, options); %New trajectory
    
    %Convergence analysis
    if (norm(error) < tol)
        GoOn = false;
    else
        iter = iter+1;
    end
end

%Final trajectory 
St = S;

dV0(1:3,1) = dV(:,iter-1);                    %Initial rendezvous impulse 
dVf(1:3,1) = -S(end,10:12).';                 %Final rendezvous impulse 

%Total maneuver metrics 
dV1(1:3,1) = dV0(:,1)+dVf(:,1);               %L1 norm of the impulses 
dV2(1) = norm(dV0(:,1))+norm(dVf(:,1));       %L2 norm of the impulses 

Pass(1) = ~GoOn;

%% GNC: one impulsive rendezvous, single shooting scheme %% 
%Differential corrector set up
S = Sn(1:index,:);                  %Reinitiate the trajectory
s0 = Sn(1,:).';                     %Reinitiate initial conditions
maxIter = 20;                       %Maximum number of iterations
tol = 1e-10;                        %Differential corrector tolerance
GoOn = true;                        %Convergence boolean 
iter = 1;                           %Initial iteration 

%Preallocation 
dV = zeros(3,maxIter);              %Targeting impulse

%First impulse: targeting 
while ((GoOn) && (iter < maxIter))
    %Initial impulse
    error = S(end,7:12).';                                                  %Velocity error
    STM = reshape(S(end,13:end), [length(r_t0) length(r_t0)]);              %STM evaluated at time tf
    dV(:,iter) = (STM(:,4:6).'*STM(:,4:6))^(-1)*STM(:,4:6).'*error;         %Required impulse
    
    %Dynamics integration 
    s0(10:12) = s0(10:12)-dV(:,iter);                                                       %New initial conditions
    [~, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke V', t, s), tspan, s0, options); %New trajectory
    
    %Convergence analysis
    if (norm(error) < tol)
        GoOn = false;
    else
        iter = iter+1;
    end
end

dV0(:,2) = dV(:,iter-1);                 %Initial rendezvous impulse 
dVf(:,2) = zeros(3,1);                   %Final rendezvous impulse 

%Total maneuver metrics 
dV1(1:3,1) = dV0(:,2)+dVf(:,2);          %L1 norm of the impulses 
dV2(1) = norm(dV0(:,2))+norm(dVf(:,2));  %L2 norm of the impulses 

Pass(2) = ~GoOn;

%% Results %% 
%Print results 
disp('SIMULATION RESULTS: ')
if (Pass(1))
    disp('   Two impulsive rendezvous was achieved');
    fprintf('   Initial impulse: %.4ei %.4ej %.4ek \n', dV0(1,1), dV0(2,1), dV0(3,1));
    fprintf('   Final impulse: %.4ei %.4ej %.4ek \n', dVf(1,1), dVf(2,1), dVf(3,1));
    fprintf('   Delta V budget (L1 norm): %.4ei %.4ej %.4ek \n', dV1(1,1), dV1(2,1), dV1(3,1));
    fprintf('   Delta V budget (L2 norm): %.4e \n', dV2(:,1));
else
    disp('    Two impulsive rendezvous was not achieved');
end

disp(' ');
if (Pass(2))
    disp('   One impulsive rendezvous was achieved.');
    fprintf('   Initial impulse: %.4ei %.4ej %.4ek \n', dV0(1,2), dV0(2,2), dV0(3,2));
    fprintf('   Final impulse: %.4ei %.4ej %.4ek+\n', dVf(1,2), dVf(2,2), dVf(3,2));
    fprintf('   Delta V budget (L1 norm): %.4ei %.4ej %.4ek \n', dV1(1,2), dV1(2,2), dV1(3,2));
    fprintf('   Delta V budget (L2 norm): %.4e', dV2(:,2));
else
     disp('   One impulsive rendezvous was not achieved.');
end

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
plot3(St(:,7), St(:,8), St(:,9)); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative motion in the configuration space');

%%
%Rendezvous animation 
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

