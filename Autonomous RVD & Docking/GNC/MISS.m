%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 9: Multiple impulse rendezvous using single shooting STM %% 
% This script provides an interface to test the multi-impusilve rendezvous strategy using the STM of 
% the dynamical model. 

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% The final objective is to provide a complete rendezvous condition

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
maxIter = 50;                       %Maximum number of iterations
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
r_t0 = target_orbit.Trajectory(10,1:6);                    %Initial target conditions
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

%% GNC: multi impulsive rendezvous, generalized targetting approach %%
%Differential corrector set up
maxIter = 200;                              %Maximum number of iterations
tol = 1e-5;                                 %Differential corrector tolerance
S = S(1:index,:);                           %Restrict the time integration span
T = index*dt;                               %Flight time along the arc
nodes = 2;                                  %Number of nodes to compute
GoOn = true;                                %Convergence boolean 
iter = 1;                                   %Initial iteration 

%Select impulsive times 
impulses = 5;                               %Number of impulses to be made
tspan = 0:dt:T;                             %Integration time span
if (impulses ~= 0)
    times = T*rand(1,impulses);             %Times to impulse the spacecraft
    times = fix(times/dt);                  %Position along the time span to impulse the spacecraft
    times = sort(times);                    %Sort the times at which the impulses are made
    times = [1 times size(S,1)];            %Impulses time
    impulses = length(times);               %Update the number of impulses to take into accoun the initial and final ones
else
    times = [1 size(S,1)];                  %Simply two impulsive rendezvous strategy
    impulses = length(times);               %Update the number of impulses to take into accoun the initial and final ones
end

%Cost function 
cost = 'Position';                          %Make impulses to target the state

%Implementation 
while ((GoOn) && (iter < maxIter))    
    %Recompute initial conditions
    switch (cost)
        case 'Position' 
            xf = S(end,7:9);                %Final positon state     
            STM = zeros(3,3*impulses);      %Preallocation of the STM
            
            %Compute the STM 
            for i = 1:length(times)
                aux = reshape(S(times(i),13:end), [length(r_t0) length(r_t0)]);      %STM evaluated at time ti
                STM(:,1+3*(i-1):3*i) = aux(1:3,4:6);
            end
                        
        case 'State' 
            xf = S(end,7:12);                %Final state
            STM = zeros(6,6*impulses);       %Preallocation of the STM
            
            %Compute the STM 
            for i = 1:length(times)
                STM(:,1+6*(i-1):6*i) = reshape(S(times(i),13:end), [length(r_t0) length(r_t0)]);    %STM evaluated at time ti
            end
            STM = STM(:,4:6); 
            
        otherwise
            error('No valid cost function was selected');
    end
    
    %Compute the error and the impulses 
    dV = -pinv(STM)*xf.';                           %Impulses
    
    %Reintegrate the trajectory
    for i = 1:length(times)-1
        %Integration time span
        Dt = tspan(times(i)):dt:tspan(times(i+1));
        
        %New trajectory
        [~, s] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke V', t, s), Dt, S(times(i),:), options); 
        
        %Update new initial conditions
        S(times(i):times(i+1),:) = s;
        
        %Make the impulse
        S(times(i),10:12) = S(times(i),10:12) + dV(1+3*(i-1):3*i,:).';
    end
    
    %Include the final impulse 
    S(end,10:12) = S(end,10:12) + dV(end-2:end,:).';
    
    %Convergence analysis 
    if (norm(S(end,7:12)) < tol)
        GoOn = false;                       %Stop the method
    else
        iter = iter+1;                      %Update the iterations
    end
end

%Output 
St = S; 

%Total maneuver metrics 
dV1(1:3,1) = sum(dV,1);     %L1 norm of the impulses 
dV2(1) = norm(dV);          %L2 norm of the impulses 

Pass = ~GoOn;

%% Results %% 
disp('SIMULATION RESULTS: ')
if (Pass)
    disp('   Multi impulsive rendezvous was achieved');
    fprintf('   Delta V budget (L1 norm): %.4ei %.4ej %.4ek \n', dV1(1,1), dV1(2,1), dV1(3,1));
    fprintf('   Delta V budget (L2 norm): %.4e \n', dV2(:,1));
else
    disp('    Multi impulsive rendezvous was not achieved');
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