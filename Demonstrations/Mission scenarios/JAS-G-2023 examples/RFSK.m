%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 15/01/23 % 

%% Relative Floquet Stationkeeping demonstration %% 
% This script provides an interface to test the RFSK strategies for optimal long term stationkeeping

% Units are non-dimensional and solutions are expressed in the synodic
% reference frame as defined by Howell, 1984.

%% Set up %%
% Set up graphics 
set_graphics();

% Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
% Phase space dimension 
n = 6; 

% Time span 
dt = 1e-3;                          % Time step
tf = 1.2*pi;                        % Rendezvous time

% CR3BP constants 
mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
L = libration_points(mu);           % System libration points
Lem = 384400e3;                     % Mean distance from the Earth to the Moon
T0 = 28*86400/(2*pi);               % Characteristic time of the Earth-Moon system
Vc = 1.025e3;                       % Characteristic velocity of the Earth-Moon system

% Differential corrector set up
nodes = 10;                         % Number of nodes for the multiple shooting corrector
maxIter = 20;                       % Maximum number of iterations
tol = 1e-10;                        % Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
% Halo characteristics 
Az = 30e6;                                                          % Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
Ln = 2;                                                             % Orbits around L1
gamma = L(end,Ln);                                                  % Li distance to the second primary
m = 1;                                                              % Number of periods to compute

% Compute a halo seed 
halo_param = [-1 Az Ln gamma m];                                    % Southern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          % Generate a halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Reference orbit %%
% Integration of the model
s0 = target_orbit.Trajectory(1,1:6);
s0 = [s0 reshape(eye(n), [1 n^2])];
tspan = 0:dt:target_orbit.Period;
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);
Sn = repmat(S, floor(tf/target_orbit.Period)+1, 1);     

% Guidance law for the Jacobi constant
Sg = jacobi_constant(mu, Sn(1,1:n).');

% Guidance law coefficients
order = 30;                                                 % Orbit approximation order
[Cp, Cv, Cg] = CTR_guidance(order, tspan, S(:,1:6));       % Guidance law coefficients

%% GNC algorithms definition 
GNC.Algorithms.Guidance = '';                         % Guidance algorithm
GNC.Algorithms.Navigation = '';                       % Navigation algorithm
GNC.Algorithms.Control = 'RFSK';                      % Control algorithm

GNC.Guidance.Dimension = 9;                           % Dimension of the guidance law
GNC.Control.Dimension = 3;                            % Dimension of the control law

GNC.Navigation.NoiseVariance = 0;

GNC.System.mu = mu;                                   % Systems's reduced gravitational parameter
GNC.Control.RFSK.method = 'SDRE';                     % Solver to be used
GNC.Control.RFSK.Q = 1*eye(2);                        % Penalty on the state error
GNC.Control.RFSK.M = 5e-3*eye(3);                     % Penalty on the control effort
GNC.Control.RFSK.Reference = Sg;                      % Reference state
GNC.Control.RFSK.Period = target_orbit.Period;        % Target orbital period

GNC.Tmax = 0.5e-3*(T0^2/Lem);

% Floquet mode matrix
[E, lambda] = eig(reshape(Sn(end,n+1:end), [n n]));
GNC.Control.RFSK.FloquetModes = diag(log(diag(lambda))/target_orbit.Period);
for i = 1:size(E,2)
    E(:,i) = E(:,i)/lambda(i,i);
end
GNC.Control.RFSK.FloquetDirections = E; 

% Final LQR controller
A = [lambda(1,1) 0; 0 0];                              % State matrix
dJ = jacobi_gradient(mu, Sn(1,1:n).').';               % Jacobi gradient 
C = dJ*E*lambda;                                       % Gradient of the Jacobi constant with the Floquet variables
A(2,1) = C(1);                                         % Final state matrix
B = E^(-1)*[zeros(3);eye(3)];                          % Control input matrix
B = [B(1,:); dJ*E*B];                                  % Final control input matrix

GNC.Control.RFSK.K = lqr(A,B,GNC.Control.RFSK.Q,GNC.Control.RFSK.M);                   
K = 0.3; 

% HSK controller 
GNC.Algorithms.Guidance = 'CTR';                %Guidance algorithm
GNC.Algorithms.Navigation = '';                 %Navigation algorithm
GNC.Algorithms.Control = 'HSK';                 %Control algorithm

GNC.Guidance.Dimension = 9;                     %Dimension of the guidance law
GNC.Control.Dimension = 3;                      %Dimension of the control law

GNC.Guidance.CTR.Order = order;                 %Order of the approximation
GNC.Guidance.CTR.TOF = tspan(end);              %Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;     %Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg; %Coefficients of the Chebyshev approximation

GNC.System.mu = mu;                             %Systems's reduced gravitational parameter
GNC.Control.HSK.Q = eye(1);                     %Penalty on the state error
GNC.Control.HSK.M = eye(3);                     %Penalty on the control effort

%% GNC: RFSK stationkeeping control law
% Initial conditions 
k = dimensionalizer(Lem, 1, 1, 190e3, 'Position', 0);     % Noise gain
k = [repmat(k,1,3) 4*ones(1,3)/Vc];
r_t0 = target_orbit.Trajectory(1,1:6);                    % Initial guidance target conditions

s0 = [192e3 -192e3 192e3 13 0.5 0.5]./[Lem Lem Lem Vc Vc Vc];

% Compute the trajectory
s0 = [r_t0 s0 reshape(eye(n), [1 n^2])];

% Nominal trajectory
[~, Sr] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0(1:2*n), options);

% Stationkeeping trajectory
iter = 1; 
Time = zeros(1,iter);
for i = 1:iter
    tic
    [~, Staux1] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan(tspan < K*target_orbit.Period), s0, options);                   
    [~, Staux2] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan(tspan >= K*target_orbit.Period), Staux1(end,:), options);
    Time(i) = toc;
end
St = [Staux1; Staux2];
Time = mean(Time);

Error in time 
[~, merit(:,1)] = figures_merit(tspan, St(:,n+1:2*n));
[e, merit(:,2)] = figures_merit(tspan, Sr(:,n+1:2*n));

% Control law
[~, ~, u] = GNC_handler(GNC, Staux1(:,1:n), Staux1(:,n+1:end), tspan(tspan < K*target_orbit.Period));

% Control integrals
effort = control_effort(tspan(tspan < K*target_orbit.Period), u, false)*Vc;
u_m(1) = Lem/T0^2*min(sqrt(dot(u,u,1)))*1E3;
u_m(2) = Lem/T0^2*max(sqrt(dot(u,u,1)))*1E3;

% Ratio
a = 1-norm(St(end,7:9))/norm(St(1,7:9));

% Final absolute trajectory
St = St(:,1:n)+St(:,n+1:2*n);
Sr = Sr(:,1:n)+Sr(:,n+1:2*n);

%% Results %% 
% Orbit plotting
figure(1) 
view(3) 
hold on
index = floor(linspace(1,size(St,1),50));
scatter3(Sn(index,1), Sn(index,2), Sn(index,3), 'b','filled'); 
plot3(St(:,1), St(:,2), St(:,3), 'k','LineWidth', 0.9); 
%plot3(Sr(:,1), Sr(:,2), Sr(:,3), 'r', 'LineWidth', 0.9); 
hold off
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;
legend('Reference orbit', 'Stationkeeping orbit')
% legend('Reference orbit', 'Stationkeeping orbit', 'Natural orbit')

% Phase-space evolution
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

% Error plot
figure(4)
plot(tspan, log(e)); 
xlabel('Nondimensional epoch');
ylabel('Absolute error $\log{e}$');
grid on;
title('Absolute rendezvous error in the relative space');

figure 
plot(tspan(tspan < K*target_orbit.Period), (sqrt(dot(u,u,1))*Lem/T0^2*1E3))
xlabel('$t$')
ylabel('$||\mathbf{u}||$')
grid on;

% Rendezvous animation 
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
