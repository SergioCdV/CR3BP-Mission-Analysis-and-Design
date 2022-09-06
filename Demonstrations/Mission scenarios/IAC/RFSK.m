%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 27/08/22 % 

%% Relative Floquet Stationkeeping demonstration for IAC 2022 %% 
% This script provides an interface to test the RFSK strategies for optimal long term stationkeeping.

% The relative motion of two spacecraft in the same halo orbit around L1 in the
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

%Time span 
dt = 1e-3;                          %Time step
tf = 1.5*pi;                          %Rendezvous time

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
Az = 35e6;                                                          %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Modelling in the synodic frame %%
%Integration of the model
s0 = target_orbit.Trajectory(1,1:6);
s0 = [s0 reshape(eye(n), [1 n^2])];
tspan = 0:dt:target_orbit.Period;
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);
Sn = repmat(S, floor(tf/target_orbit.Period)+1, 1);     

%Guidance law coefficients
Sg = jacobi_constant(mu,Sn(1,1:n).');

%% GNC algorithms definition 
GNC.Algorithms.Guidance = '';                         %Guidance algorithm
GNC.Algorithms.Navigation = '';                       %Navigation algorithm
GNC.Algorithms.Control = 'MSKLQR';                    %Control algorithm

GNC.Guidance.Dimension = 9;                           %Dimension of the guidance law
GNC.Control.Dimension = 3;                            %Dimension of the control law

GNC.Navigation.NoiseVariance = 0;                     %Noise variance

GNC.System.mu = mu;                                   %Systems's reduced gravitational parameter
GNC.Control.MSKDRE.Q = 1e3*eye(2);                    %Penalty on the state error
GNC.Control.MSKDRE.M = eye(3);                        %Penalty on the control effort
GNC.Control.MSKDRE.Reference = Sg;                    %Reference state
GNC.Control.MSKDRE.Period = target_orbit.Period;      %Target orbital period

%Floquet mode matrix
[E, lambda] = eig(reshape(Sn(end,n+1:end), [n n]));
GNC.Control.MSKDRE.FloquetModes = diag(log(diag(lambda))/target_orbit.Period);
for i = 1:size(E,2)
    E(:,i) = E(:,i)/lambda(i,i);
end
GNC.Control.MSKDRE.FloquetDirections = E; 

%Constant controller
A = [lambda(1,1) 0; 0 0];                              %State matrix
dJ = jacobi_gradient(mu, Sn(1,1:n).').';               %Jacobi gradient 
C = dJ*E*lambda;                                       %Gradient of the Jacobi constant with the Floquet variables
A(2,1) = C(1);                                         %Final state matrix
B = E^(-1)*[zeros(3);eye(3)];                          %Control input matrix
B = [B(1,:); dJ*B];                                    %Final control input matrix

%Final LQR controller
GNC.Control.MSKLQR.K = lqr(A,B,GNC.Control.MSKDRE.Q,GNC.Control.MSKDRE.M);                   

%% GNC: MLQR control law
%Initial conditions 
k = dimensionalizer(Lem, 1, 1, 5e3, 'Position', 0);         %Noise gain
k = [repmat(k,1,3) repmat(k,1,3)/1e6];
r_t0 = target_orbit.Trajectory(1,1:6);                      %Initial guidance target conditions
s0 = r_t0+k*normrnd(0,1);                                   %Noisy initial conditions
tspan = 0:dt:tf;                                            %New integration time span

rev = floor(tf/target_orbit.Period);

%Compute the trajectory
s0 = [r_t0 s0-r_t0 reshape(eye(n), [1 n^2])];
[~, Sr] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
tic
[~, Staux1] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan(tspan < 0.3), s0, options);
[~, Staux2] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan(size(Staux1,1):end), Staux1(end,:), options);
toc
St = [Staux1; Staux2(2:end,:)];

%Error in time 
[e, merit] = figures_merit(tspan, [St(:,1:n) abs(St(:,1:n)-Sn(1:size(St,1),1:n))]);

%Control law
[~, ~, u] = GNCc_handler(GNC, Staux1(:,1:n), Staux1(:,n+1:end), tspan(tspan < 0.3));

%Control integrals
energy = control_effort(tspan(tspan < 0.3), u, false);

%Unstable component
alpha = zeros(1,size(St,1)); 
for i = 1:size(St,1)
    STM = reshape(St(i,13:end), [n n]);
    [V,lambda] = eig(STM); 
    Alpha = V^(-1)*St(i,7:12).'; 
    alpha(i) = Alpha(1);
end

%Final absolute trajectory
St = St(:,1:n)+St(:,n+1:2*n);
Sr = Sr(:,1:n)+Sr(:,n+1:2*n);

%% Results %% 
%Plot results 
figure(1) 
view(3) 
hold on
index = floor(linspace(1,size(St,1),50));
scatter3(Sn(index,1), Sn(index,2), Sn(index,3), 'b','filled'); 
plot3(Sr(:,1), Sr(:,2), Sr(:,3), 'r', 'LineWidth', 0.9); 
plot3(St(:,1), St(:,2), St(:,3), 'k','LineWidth', 0.9); 
hold off
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;
legend('Reference orbit', 'Natural orbit', 'Stationkeeping orbit')
axis('equal')

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan, St(:,1)); 
plot(tspan, St(:,2)); 
plot(tspan, St(:,3)); 
hold off
xlabel('$t$');
ylabel('$\mathbf{\rho}$');
grid on;
legend('$x$', '$y$', '$z$');
title('Position in time');
subplot(1,2,2)
hold on
plot(tspan, St(:,4)); 
plot(tspan, St(:,5)); 
plot(tspan, St(:,6)); 
hold off
xlabel('$t$');
ylabel('$\dot{\mathbf{\rho}}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');

%Error plot
figure(4)
plot(tspan, log(e)); 
xlabel('Nondimensional epoch');
ylabel('Absolute error $\log{e}$');
grid on;
title('Absolute rendezvous error in the relative space');

figure 
plot(tspan(tspan < 0.1), abs(alpha(tspan < 0.1)))
xlabel('$t$')
ylabel('$\alpha_1$')
grid on; 

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
        drawnow;
        delete(T); 
        delete(V);
    end
    hold off
end