%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 12: APF guidance-control law %% 
% This script provides an interface to test Artificial Potential Functions guidance laws for
% rendezvous missions.

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
Az = 200e6;                                                         %Orbit amplitude out of the synodic plane. 
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
num = 2;                                                            %Number of orbits to continuate
method = 'SPC';                                                     %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                        %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};                 %Object and characteristics to continuate
corrector = 'Plane Symmetric';                                      %Differential corrector method
direction = 1;                                                      %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                                 %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(2,:), maxIter, tol);

%% Modelling in the synodic frame %%
index = fix(tf/dt);                                         %Rendezvous point
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% APF guidance scheme
%Compute some random objects  in the relative phase space 
So = rand(3,3);

%Compute the guidance law
Sg = apf_guidance('Steady', false, Sn(1,7:12), So, dt, zeros(3,1));

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

%Configuration space error 
figure(3)
plot(tspan, e); 
xlabel('Nondimensional epoch');
ylabel('Absolute error');
grid on;
title('Absolute error in the configuration space (L2 norm)');

%Rendezvous animation 
if (false)
    figure(4) 
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

%% Auxiliary functions
%Guidance law function 
function [Sg, phi] = apf_guidance(dynamics, safe_corridor, state, obstacles_states, dt, Sg0)
    %Constants 
    m = 6;              %State dimension 
    Q = eye(m/2);       %Penalty on the distance to the origin
    R = eye(m/2);       %Penalty on the distance to the obstacles 
    
    chi = deg2rad(45);  %Safety corridor angle
    rho = 1e-4;         %Safety distance to the docking port
    
    %Integration options 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);   
    
    %Sanity check on dimension 
    if (size(state,1) ~= m)
        state = state.';
    end
    
    %Compute the attractive APF value
    switch (dynamics)
        case 'Steady'
            phi = (1/2) * state(1:3).'*Q*state(1:3);        %Attractive steady APF
            dPhi = Q*state(1:3);                            %Gradient of the APF 
            hPhi = Q;                                       %Hessian of the APF 
        case 'Unsteady'
        otherwise
            error('No valid APF dynamics were selected');
    end

    %Compute the repulsive APF values
    for j = 1:size(obstacles_states,2)
        %Relative distance to each obstacle
        rel_state = state(1:3)-obstacles_states(1:3,j);       

        %Repulsive APF
        phi = phi + (1/2) * (state(1:3).'*Q*state(1:3))/(rel_state.'*R*rel_state-1); 
        dPhi = dPhi + (rel_state.'*R*rel_state-1)^(-2)*((rel_state.'*R*rel_state-1)*Q*state(1:3) ...
                                                        -(state(1:3).'*Q*state(1:3))*R*rel_state);
        hPhi = hPhi - 2*((rel_state.'*R*rel_state-1)*Q*state(1:3)-(state(1:3).'*Q*state(1:3))*R*rel_state)*...
                       (rel_state.'*R*rel_state-1)^(-3*)*R*rel_state
    end

    %Compute the safety APF value
    if (safe_corridor)
        f = state(2)^2+state(3)^2-((state(1)^2*tan(chi))/(2*rho-state(1)));     %Corridor surface
        phi = phi + exp(-(1/2)*norm(f));                                        %Safety APF
        dPhi = dPhi + (1/2) * exp(-(1/2)*norm(f)); 
    end

    %Compute the guidance trajectory
    [~,r] = ode45(@(t,r)(-dPhi), [0 dt], Sg0, options);         %Reference position integration 
    Sg(1:3) = r(end,:).';                                       %Reference position
    Sg(4:6) = -dPhi;                                            %Reference velocity
    Sg(7:9) = -dPhi.'*hPhi;                                     %Reference acceleration
    
    %Sanity check on dimension 
    Sg = Sg.'; 
end