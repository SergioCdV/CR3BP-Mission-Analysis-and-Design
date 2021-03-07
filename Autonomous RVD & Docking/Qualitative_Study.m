%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 20/02/21 % 

%% First contact %% 
% This scripts provides a first contact with the RVD project. 

% The relative motion of two spacecraft in a halo orbit around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion % 

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode45)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
% Time span 
dt = 1e-3;                  %Time step
tmax = 2*pi;                %Maximum time of integration (corresponding to a synodic period)
tspan = 0:dt:tmax;          %Integration time span

% CR3BP constants 
mu = 0.0121505;             %Earth-Moon reduced gravitational parameter
L = libration_points(mu);   %System libration points
Lem = 384400e3;             %Mean distance from the Earth to the Moon

% CR3BP integration flags 
flagVar = 1;                %Integrate the dynamics and the first variotional equations 
direction = 1;              %Integrate forward in time 

% Differential corrector set up
nodes = 10;                 %Number of nodes for the multiple shooting corrector
maxIter = 20;               %Maximum number of iterations
tol = 1e-7;                 %Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
%Halo characteristics 
Az = 75e6;                                                  %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
K = 1e-7;                                                   %Nondimensional distance from target to chaser

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[halo_orbit, state] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Modelling %% 
r_t0 = halo_orbit.Trajectory(1,1:6);        %Initial target conditions
r_c0(1:6) = r_t0(1:6)+K*rand(1,6);          %Initial chaser conditions 
rho0 = K*rand(1,6);                         %Initial relative conditions
s0 = [r_t0 rho0];                           %Initial conditions of the target and the relative state

%Integration of the model
[~, S_c] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, r_c0, options);
[t, S] = ode113(@(t,s)full_model(mu, true, false, t, s), tspan, s0, options);

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                  %Reconstructed chaser motion

%Compute the error 
error = S_c-S_rc;

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

figure(2) 
plot(t, error); 
xlabel('Nondimensional epoch'); 
ylabel('Relative error'); 

%% Relative motion model 
function [drho] = relative_motion(mu, s_t, s_r)
    %Constants of the system 
    mu1 = 1-mu;             %Reduced gravitational parameter of the first primary 
    mu2 = mu;               %Reduced gravitational parameter of the second primary 
    
    %State variables 
    r_t = s_t(1:3);         %Synodic position of the target
    r_r = s_r(1:3);         %Synodic relative position 
    v_r = s_r(4:6);         %Synodic relative velotice 
    
    %Synodic position of the primaries 
    R1 = [-mu; 0; 0];       %Synodic position of the first primary
    R2 = [1-mu; 0; 0];      %Synodic position of the second primary
    
    %Equations of motion 
    drho = [v_r; 
            r_r(1)+2*v_r(2)-mu1*((r_t(1)+mu)/norm(r_t-R1)^3-(r_t(1)+r_r(1)+mu)/norm(r_t-R1+r_r)^3)-mu2*((r_t(1)-(1-mu))/norm(r_t-R2)^3-(r_t(1)+r_r(1)-(1-mu))/norm(r_t-R2+r_r)^3); 
            r_r(2)-2*v_r(1)-mu1*(r_t(2)/norm(r_t-R1)^3-(r_t(2)+r_r(2))/norm(r_t-R1+r_r)^3)-mu2*(r_t(2)/norm(r_t-R2)^3-(r_t(2)+r_r(2))/norm(r_t-R2+r_r)^3); 
            -mu1*(r_t(3)/norm(r_t-R1)^3-(r_t(3)+r_r(3))/norm(r_t-R1+r_r)^3)-mu2*(r_t(3)/norm(r_t-R2)^3-(r_t(3)+r_r(3))/norm(r_t-R2+r_r)^3)];
end

function [ds] = full_model(mu, true, false, t, s)
    %State variables 
    s_t = s(1:6);       %State of the target
    rho = s(7:12);      %State of the chaser
    
    %Equations of motion 
    ds_t = cr3bp_equations(mu, true, false, t, s_t);        %Target equations of motion
    drho = relative_motion(mu, s_t, rho);                   %Relative motion equations
    
    %Vector field 
    ds = [ds_t; drho];
end