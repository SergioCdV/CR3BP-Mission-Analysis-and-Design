%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 9: Multiple impulse rendezvous using the SMT, single shooting %% 
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
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state
Phi = eye(length(r_t0));                                    %Initial STM
Phi = reshape(Phi, [length(r_t0)^2 1]);                     %Initial STM
s0 = [s0; Phi];                                             %Initial conditions

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspann, s0, options);
Sn = S;              

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% GNC: multi impulsive rendezvous, generalized targetting approach %%
%Differential corrector set up
maxIter = 200;                              %Maximum number of iterations
tol = 1e-8;                                 %Differential corrector tolerance
S = S(1:index,:);                           %Restrict the time integration span
T = index*dt;                               %Flight time along the arc
nodes = 2;                                  %Number of nodes to compute
GoOn = true;                                %Convergence boolean 
iter = 1;                                   %Initial iteration 

%Select impulsive times 
impulses = 3;                               %Number of impulses to be made
if (impulses ~= 0)
    times = T*rand(1,impulses);             %Times to impulse the spacecraft
    times = fix(times/dt);                  %Position along the time span to impulse the spacecraft
    times = sort(times);                    %Sort the times at which the impulses are made
    times = [1 times size(S,1)-1];          %Impulses time
    impulses = length(times);               %Update the number of impulses to take into accoun the initial and final ones
else
    times = [2 size(S,1)-1];                %Simply two impulsive rendezvous strategy
    impulses = length(times);               %Update the number of impulses to take into accoun the initial and final ones
end

%Cost function 
cost = 'Position';                          %Make impulses to target the state

%Initial conditions 
S0 = s0;

%Preallocation of the impulses 
dV = zeros(3*impulses, maxIter);

%Implementation 
while ((GoOn) && (iter < maxIter))    
    %Recompute initial conditions
    switch (cost)
        case 'Position' 
            xf = S(end,7:9);                %Final positon state     
            STM = zeros(3,3*impulses);      %Preallocation of the STM
            
            %Compute the STM 
            for i = 1:length(times)
                STMf = reshape(S(end,13:end), [length(r_t0) length(r_t0)]);         %STM evaluated at time tf
                STM0 = reshape(S(times(i),13:end), [length(r_t0) length(r_t0)]);    %STM evaluated at time ti
                STMdt = STMf*STM0^(-1);                                             %STM evaluated between times tf-ti
                STM(:,1+3*(i-1):3*i) = STMdt(1:3,4:6);
            end

        case 'Velocity' 
            xf = S(end,10:12);               %Final state
            STM = zeros(3,3*impulses);       %Preallocation of the STM
            
            %Compute the STM 
            for i = 1:length(times)
                STMf = reshape(S(end,13:end), [length(r_t0) length(r_t0)]);         %STM evaluated at time tf
                STM0 = reshape(S(times(i),13:end), [length(r_t0) length(r_t0)]);    %STM evaluated at time ti
                STMdt = STMf*STM0^(-1);                                             %STM evaluated between times tf-ti
                STM(:,1+3*(i-1):3*i) = STMdt(4:6,4:6);
            end
                        
        case 'State' 
            xf = S(end,7:12);                %Final state
            STM = zeros(6,3*impulses);       %Preallocation of the STM
            
            %Compute the STM 
            for i = 1:length(times)
                STMf = reshape(S(end,13:end), [length(r_t0) length(r_t0)]);         %STM evaluated at time tf
                STM0 = reshape(S(times(i),13:end), [length(r_t0) length(r_t0)]);    %STM evaluated at time ti
                STMdt = STMf*STM0^(-1);                                             %STM evaluated between times tf-ti
                STM(:,1+3*(i-1):3*i) = STMdt(:,4:6);
            end
                        
        otherwise
            error('No valid cost function was selected');
    end
    
    %Compute the error and the impulses 
    if (size(STM,1) == length(xf))
        dV(:,iter) = -STM\xf.';                      %Impulses
    else
        dV(:,iter) = -STM.'*(STM*STM.')^(-1)*xf.';   %Impulses
    end
    
    %Reintegrate the trajectory
    for i = 1:length(times)               
        %Make the impulse
        S(times(i),10:12) = S(times(i),10:12) + dV(1+3*(i-1):3*i,iter).';
        
        %Integration time span
        if (i ~= length(times))
            Dt = tspan(times(i)):dt:tspan(times(i+1));
        else
            Dt = tspan(times(i)):dt:tspan(end);
        end
        
        %New trajectory
        [~, s] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), Dt, S(times(i),:), options); 
        
        %Update new initial conditions
        S(times(i):times(i)+size(s,1)-1,:) = s;
    end
        
    %Convergence analysis 
    if (norm(xf) < tol)
        GoOn = false;                       %Stop the method
    else
        iter = iter+1;                      %Update the iterations
    end
end

%Output 
St = S; 
dV(end+1:end+3,iter) = -St(end,10:12).';  %Docking impulse
St(end,10:12) = zeros(1,3);               %Null relative velocity state

%Total maneuver metrics 
dV = sum(dV,2);
dV = reshape(dV, [3 impulses+1]);         %Reshape the impulses array
dV1(1:3,1) = sum(dV,2);                   %L1 norm of the impulses 
dV2(1) = sum(sqrt(dot(dV,dV,2)));         %L2 norm of the impulses 

Pass = ~GoOn;

%Error in time 
e = zeros(1,size(St,1));                  %Preallocation of the error
for i = 1:size(St,1)
    e(i) = norm(St(i,7:12));
end
e(1) = norm(S0(7:12));                    %Initial error before the burn

%Compute the error figures of merit 
ISE = trapz(tspan, e.^2);
IAE = trapz(tspan, abs(e));

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

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan, St(:,7)); 
plot(tspan, St(:,8)); 
plot(tspan, St(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('x coordinate', 'y coordinate', 'z coordinate');
title('Relative position evolution');
subplot(1,2,2)
hold on
plot(tspan, St(:,10)); 
plot(tspan, St(:,11)); 
plot(tspan, St(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinate');
grid on;
legend('x velocity', 'y velocity', 'z velocity');
title('Relative velocity evolution');

%Configuration space error 
figure(4)
plot(tspan, e); 
xlabel('Nondimensional epoch');
ylabel('Absolute error');
grid on;
title('Absolute error in the configuration space (L2 norm)');

%Rendezvous animation 
if (false)
    figure(5) 
    view(3) 
    grid on;
    hold on
    plot3(St(1:index,1), St(1:index,2), St(1:index,3), 'k-.'); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    title('Rendezvous simulation');
    for i = 1:size(Sc,1)
        T = scatter3(St(i,1), St(i,2), St(i,3), 30, 'b'); 
        V = scatter3(St(i,1)+St(i,7), St(i,2)+St(i,8), St(i,3)+St(i,9), 30, 'r');

        drawnow;
        delete(T); 
        delete(V);
    end
    hold off
end