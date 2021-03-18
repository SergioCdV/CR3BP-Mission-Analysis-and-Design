%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 20/02/21 % 

%% Linear relative motion models %% 
% This scripts provides a complete script dedicated to linear models 
% for relative motion in the CR3BP.

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
dt = 1e-3;                          %Time step
tmax = 0.8*pi;                        %Maximum time of integration (corresponding to a synodic period)
tspan = 0:dt:tmax;                  %Integration time span

% CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
L = libration_points(mu);           %System libration points
Lem = 384400e3;                     %Mean distance from the Earth to the Moon

% CR3BP integration flags 
flagVar = 1;                        %Integrate the dynamics and the first variotional equations 
direction = 1;                      %Integrate forward in time 

% Differential corrector set up
nodes = 10;                         %Number of nodes for the multiple shooting corrector
maxIter = 20;                       %Maximum number of iterations
tol = 1e-10;                        %Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
%Halo characteristics 
Az = 200e6;                                                 %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
cn = legendre_coefficients(mu, Ln, gamma, 2);               %Legendre coefficients of interest                                   
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
[~, S_cf] = ode113(@(t,s)fulrel_motion(mu, true, false, true, t, s), tspan, s0, options);
[t, S] = ode113(@(t,s)linrel_model(mu, true, false, 'Libration', t, s, cn(2)), tspan, s0, options);

%Reconstructed chaser motion 
S_cf = S_cf(:,1:6)+S_cf(:,7:12);                            %Reconstructed chaser motion via linear models
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via linear models
error_rel = S_cf-S_rc;                                      %Error via the Encke method and the linear models
error_abs = S_c-S_rc;                                       %Error between the full integration and the linear models

er = zeros(size(error_rel,1), 1);                           %Position error (L2 norm) via linear models
ea = zeros(size(error_abs,1), 1);                           %Position error (L2 norm) via direct integration
for i = 1:size(error_rel,1)
    er(i) = norm(error_rel(i,1:3));
    ea(i) = norm(error_abs(i,1:3));
end

%% Comparison with the analytical solution 
%Analytical solution
Ax = -S(1,7);                                       %Relative orbit amplitude
omega = sqrt(cn(2)+2+sqrt(9*cn(2)^2-8*cn(2)))/2;    %Orbital frequency
k = (omega^2+(1-2*cn(2)))/(2*omega);                %Amplitude constraint gain
x = -Ax*cos(omega*tspan);                           %Analytical synodic x solution
y = k*Ax*sin(omega*tspan);                          %Analytical synodic y solution

z = (sqrt(cn(2))/(2*cn(2)))*[S(1,12)+S(1,9)*sqrt(cn(2)) -S(1,12)+S(1,9)*sqrt(cn(2))]*[exp(sqrt(cn(2))*tspan); ...
           exp(-sqrt(cn(2))*tspan)]; 
       
rho_a = [x.' y.' z.'];
S_a(:,1:3) = S(:,1:3)+rho_a;                    %Reconstructed chaser motion via linear models

%% Results %% 
%Plot results 
figure(1) 
view(3) 
hold on
plot3(S(:,1), S(:,2), S(:,3),'r-.'); 
plot3(S_c(:,1), S_c(:,2), S_c(:,3)); 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
plot3(S_a(:,1), S_a(:,2), S_a(:,3)); 
hold off
legend('Target motion', 'Chaser motion', 'Linear integration', 'Analytical solution'); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Reconstruction of the chaser motion');

figure(2) 
hold on
plot(t, er, 'b'); 
plot(t, ea, 'r');
hold off
grid on
xlabel('Nondimensional epoch'); 
ylabel('Relative error in the synodic frame');
title('Error in the chaser velocity (L2 norm)')
legend('Linear method w.r.t Encke', 'Linear method w.r.t integration');



