%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 25/04/21 % 

%% GNC 11: Complete rendezvous mission example 1 %% 
% This script provides an interface to test the general control scheme for a rendezvous, docking and undocking mission. 

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

%Spacecraft mass 
mass = 1e-10;

%Time span 
dt = 1e-3;                          %Time step
tf = [0.6 pi pi+0.1];               %Switching times
tspan = 0:dt:tf(3);                 %Integration time span

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
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state
Phi = eye(length(r_t0));                                    %Initial STM
Phi = reshape(Phi, [length(r_t0)^2 1]);                     %Initial STM
s0 = [s0; Phi];                                             %Initial conditions

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% First phase: long-range rendezvous using the TITA approach
%Differential corrector set up
index = fix(tf(1)/dt);                      %First event: end of the long-range rendezvous phase
tspan = 0:dt:tf(1);                         %Initial integration time
maxIter = 100;                              %Maximum number of iterations
tol = 1e-5;                                 %Differential corrector tolerance
S = S(1:index,:);                           %Restrict the time integration span
T = index*dt;                               %Flight time along the arc
GoOn = true;                                %Convergence boolean 
iter = 1;                                   %Initial iteration 

%Preallocation 
dV = zeros(3,maxIter);                      %Targeting impulse

%Cost function matrices
R = eye(3);                                 %Penalty on the impulse
Qt = eye(6);                                %Penalty on the state error
M = 0.1*eye(6);                             %Penalty on the state noise
Omegat = [zeros(3,3); eye(3)];              %Derivative of the state vector with respect to the impulse
Rr = eye(3);                                %Rotational misalignment of the thrusters

%Select measuring times 
noise = true;                               %Boolean to account for state noise
measurements = 3;                           %Number of noise measurements
times = T*rand(1,measurements);             %Times to measure the state noise
times = fix(times/dt);                      %Position along the time span to measure the state noise
times = sort(times);                        %Sort the times at which the noise measurements are taken
ns = 1e-6*ones(6,1);                        %Initial state noise 
sigma = 0.01;                               %Velocity noise dependance on the velocity impulse

%Cost function 
cost = 'Position';                          %Make impulses to target position

%Initial conditions 
S0 = s0;

%Implementation 
while ((GoOn) && (iter < maxIter))
    %Compute the complete STM
    STM = reshape(S(end,13:end), [length(r_t0) length(r_t0)]);      %STM evaluated at time tf
    
    %Propagate the error 
    if (noise)
        nSTM = zeros(6,3); 
        for i = 1:measurements 
             dumbSTM = reshape(S(times(i),13:end), [length(r_t0) length(r_t0)]);     %Noise state transition matrix
             nSTM = nSTM + dumbSTM.'*M*dumbSTM*[zeros(3,3); Rr];                     %Accumulated state transition matrix
        end
        nSTM = sigma^2*[zeros(6,3) [zeros(3,3); Rr]]*nSTM;                           %Accumulated state transition matrix
        nState = sigma*ns.'*nSTM;                                                    %Accumulated noise vector
    end

    %Recompute initial conditions
    switch (cost)
        case 'Position' 
            xf = S(end,7:9);                %Final positon state
            Phi = STM(1:3,4:6);             %Needed state transition matrix
            Q = Qt(1:3,1:3);                %Penalty on the state error 
            Omega = Omegat(4:6,:);          %Derivative of the state vector with respect to the impulse
            
            %Compute the STM 
            L = Omega.'*Phi.'*Q*Phi*Omega;  %Penalty on the state error 
            STM = R+L;                      %Jacobian of the constraint function
                        
        case 'State' 
            xf = S(end,7:12);               %Final state
            Phi = STM;                      %Needed state transition matrix
            Q = Qt;                         %Penalty on the state error
            Omega = Omegat;                 %Derivative of the state vector with respect to the impulse
            
            %Compute the STM
            L = Omega.'*Phi.'*Q*Phi*Omega;  %Penalty on the state error 
            STM = R+L;                      %Jacobian of the constraint function            
            
        otherwise
            error('No valid cost function was selected');
    end
        
    %Add some noise 
    if (noise)
        STM = STM + nSTM(4:6,1:3);                             %Noise state matrix
        error = xf*Q*Phi*Omega + nState;                       %Error state (deviation from the rendezvous condition)
        dV(:,iter) = pinv(STM)*error.';                        %Needed impulse
        s0(10:12) = s0(10:12)-dV(:,iter)+sigma*Rr*dV(:,iter);  %New initial conditions
        s0(7:12) = s0(7:12)+ns;                                %New noisy initial conditions
    else
        error = xf*Q*Phi*Omega;                                %Error state (deviation from the rendezvous condition)
        dV(:,iter) = STM\error.';                              %Needed impulse
        s0(10:12) = s0(10:12)-dV(:,iter);                      %New initial conditions 
    end
    
    %Reintegrate the trajectory
    [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); %New trajectory
    
    %Convergence analysis 
    if (norm(error) < tol)
        GoOn = false;                        %Stop the method
    else
        iter = iter+1;                       %Update the iterations
    end
end

%Output 
St1 = S;                                     %Integrated trajectory
dV0(1:3,1) = sum(dV,2);                      %Initial rendezvous impulse 
dVf(1:3,1) = -S(end,10:12).';                %Final rendezvous impulse 
St1(end,10:12) = St1(end,10:12)+dVf.';       %Docking burn condition

%% Second phase: docking and coordinated flight 
%Integration time 
tspan = 0:dt:tf(2)-tf(1); 

%Initial conditions 
s0 = St1(end,:); 

%Re-integrate trajectory
[~, St2] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke SMC', t, s, false), tspan, s0, options);

%% Third phase: undocking 
%Integration time 
tspan = 0:dt:tf(3)-tf(2)-dt; 

%Initial conditions 
s0 = St2(end,:); 

%Select the restriction level of the CAM 
lambda(1) = 1e-1;                                           %Safety distance
lambda(2) = 1e-1;                                           %Safety distance

%Compute the Floquet modes at each time instant 
Monodromy = reshape(s0(13:end), [6 6]);                     %State transition matrix at each instant        
[E, sigma] = eig(Monodromy);                                %Eigenspectrum of the STM 

%Compute the maneuver
STM = [E(:,1) E(:,3:end) -[zeros(3,3); eye(3)]];             %Linear application
error = lambda(1)*rand(6,1);                                 %State error vector
maneuver = STM.'*(STM*STM.')^(-1)*error;                     %Needed maneuver
dV = real(maneuver(end-2:end));                              %Needed change in velocity

%Integrate the CAM trajectory
s0(10:12) = s0(10:12)+dV.';                                  %Update initial conditions with the velocity change

[~, St3] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 

%% Final results 
St = [St1; St2; St3];   %Complete trajectory 
tspan = 0:dt:tf(3)+dt;  %Total integration time

%% Plotting
figure(2) 
view(3) 
hold on 
plot3(St(:,7), St(:,8), St(:,9), 'r'); 
hold off
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative motion trajectory in the configuration space');

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan(1:size(St,1)), St(:,7)); 
plot(tspan(1:size(St,1)), St(:,8)); 
plot(tspan(1:size(St,1)), St(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('x coordinate', 'y coordinate', 'z coordinate');
title('Relative position evolution');
subplot(1,2,2)
hold on
plot(tspan(1:size(St,1)), St(:,10)); 
plot(tspan(1:size(St,1)), St(:,11)); 
plot(tspan(1:size(St,1)), St(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinate');
grid on;
legend('x velocity', 'y velocity', 'z velocity');
title('Relative velocity evolution');
