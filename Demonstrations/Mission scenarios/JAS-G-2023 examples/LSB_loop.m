%% Autonomous RVD and docking in the CR3BP  %%
% Date: 06/01/23

%% Set up
set_graphics(); 
clear;
close all

%% Trajectory generation 
% CR3BP constants 
mu = 0.0121505;                         % Earth-Moon reduced gravitational parameter
L = libration_points(mu);               % System libration points
Lem = 384400e3;                         % Mean distance from the Earth to the Moon
T0 = 28*86400/(2*pi);                   % Characteristic time of the Earth-Moon system
Vc = 1.025e3;                           % Characteristic velocity of the Earth-Moon system
n = 6;                                  % Phase-space dimension

% Halo characteristics 
Az = 5e6:100e4:15e6;                    % Orbit amplitude out of the synodic plane  
Az = 5e6:100e3:35e6;                    % Orbit amplitude out of the synodic plane 
Ln = 1;                                 % Libration point index
northern_flag = 1;                      % Northern/southern familiy flag

% Differential corrector set up
maxIter = 20;                           % Maximum number of iterations
tol = 1e-10;                            % Differential corrector tolerance

gamma = L(end,Ln);                      % Li distance to the second primary
m = 1;                                  % Number of periods to compute

T = 0.5e-3*(T0^2/Lem);                  % Maximum thrust in non-dimensional units
GNC_LSB.LQR.StateMatrix = 1*eye(2);     % State error weight matrix
GNC_LSB.LQR.ControlMatrix = 1e2*eye(1); % Control effort weight matrix
GNC_LSB.BSK.K1 = 10*eye(2);                % First regulation matrix
GNC_LSB.BSK.K2 = eye(1);                % Second regulation matrix
GNC_LSB.Tmax = T;                       % Constrained acceleration
GNC_LSB.TOF = 2*pi;                     % Maneuver time
GNC_LSB.Polynomial = 'Bernstein';       % Polynomial family to be used

% Target orbit
az = dimensionalizer(Lem, 1, 1, Az(1), 'Position', 0);              % Normalize distances for the E-M system
halo_param = [northern_flag az Ln gamma m];                         % Halo parameters
[halo_seed, ~] = object_seed(mu, halo_param, 'Halo');               % Generate the halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

target_state = target_orbit.Trajectory(100,1:6);                    % Initial target conditions

options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);

%% Close-loop performance example
% Compute the halo orbit
az = Az(2)/Lem;                                                   % Normalize distances for the E-M system
az = Az(20)/Lem;                                                    % Normalize distances for the E-M system
halo_param = [northern_flag az Ln gamma m];                         % Halo parameters
[halo_seed, ~] = object_seed(mu, halo_param, 'Halo');               % Generate a halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
 
initial_state = chaser_orbit.Trajectory(1,1:6);                     % Initial chaser conditions 

iter = 1;
Time = zeros(3,iter);

% Relative guidance solution 
GNC_LSB.TOF = 1;
GNC_LSB.TOF = 2;
method = "Minimum energy";                                            % Direct prescritpion 
rho = initial_state-target_state;

for i = 1:iter
    tic;
    [Sr, u_ol, tf, ~] = LSB_guidance(mu, Ln, gamma, rho, method, GNC_LSB.TOF, GNC_LSB);  
    Time(1,i) = toc;
end

tau = linspace(0,tf,size(Sr,2));
[~, Sc] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tau, target_state, options);
C = Sc.' + Sr(1:n,:);

% Total transfer metrics 
effort(:,1) = control_effort(tau, u_ol, false) * Vc;
[~, merit(:,1)] = figures_merit(tau, Sr(1:n,:).');
        
% Minimum and maximum control 
u_excl(1) = (Lem/T0^2)*min(sqrt(dot(u_ol,u_ol,1)))*1e3;
u_excl(2) = (Lem/T0^2)*max(sqrt(dot(u_ol,u_ol,1)))*1e3;

% Definition of the CL controller
GNC.Guidance.Dimension = 9;                 % Dimension of the guidance law
GNC.Control.Dimension = 3;                  % Dimension of the control law
GNC.Navigation.NoiseVariance = 0;           % Noise variance
GNC.Algorithms.Navigation = '';             % Navigation algorithm

GNC.System.mu = mu;                         % Systems's reduced gravitational parameter
GNC.System.Libration = [Ln gamma];          % Libration point ID
GNC.Tmax = T;                               % Maximum available acceleration

% Definition of the CL controller
GNC.Algorithms.Guidance = '';               % Guidance algorithm
GNC.Algorithms.Control = 'LSB';             % Control algorithm
GNC.LSB.Parameters = GNC_LSB;               % Definition of the control law
GNC.LSB.Method = method;

options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-14);

% Close loop control 
for i = 1:iter
    tic;
    [~, Sc_mpc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tau, [target_state rho], options);
    [~, ~, u_mpc] = GNC_handler(GNC, Sc_mpc(:,1:n), Sc_mpc(:,7:end), tau);
    Time(2,i) = toc;
end

effort(:,2) = control_effort(tau, u_mpc, false) * Vc;
[~, merit(:,2)] = figures_merit(tau, Sc_mpc(:,7:12));
cum_effort = Vc * cumtrapz(tau, sqrt(dot(u_mpc,u_mpc,1)));

u_excl(3) = (Lem/T0^2)*min(sqrt(dot(u_mpc,u_mpc,1)))*1e3;
u_excl(4) = (Lem/T0^2)*max(sqrt(dot(u_mpc,u_mpc,1)))*1e3;

dV_cl(1) = effort(1,2)/effort(1,1)-1;

% Guidance 
order = 30; 
[Cp, Cv, Cg, Ci] = CTR_guidance(order, tau, Sr.');

GNC.Algorithms.Guidance = 'CTR';                             % Guidance algorithm
GNC.Guidance.CTR.Order = order;                              % Order of the approximation
GNC.Guidance.CTR.TOF = tau(end);                             % Time of flight
GNC.Guidance.CTR.PositionCoefficients = Cp;     	         % Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.VelocityCoefficients = Cv;                  % Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.AccelerationCoefficients = Cg;              % Coefficients of the Chebyshev approximation
GNC.Guidance.CTR.IntegralCoefficients = Ci;                  % Coefficients of the Chebyshev approximation

% Definition of the SMC controller
GNC.Algorithms.Control = 'SMC';                              % Control algorithm

GNC.Control.SMC.Parameters = [3*1.000000000000000 5*0.985999332287318 0.006010671478548 0.013227007322678];
GNC.Control.SMC.Parameters = [3*1.000000000000000 0.985999332287318 0.006010671478548 0.013227007322678];
                
% Close loop control 
for i = 1:iter
    tic;
    [~, Sc_smc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tau, [target_state rho], options);
    [~, ~, u_smc] = GNC_handler(GNC, Sc_smc(:,1:n), Sc_smc(:,7:end), tau);
    Time(3,i) = toc;
end

% Total transfer metrics 
effort(:,3) = control_effort(tau, u_smc, false) * Vc;
[~, merit(:,3)] = figures_merit(tau, Sc_smc(:,7:12));
cum_effort(2,:) = Vc * cumtrapz(tau, sqrt(dot(u_smc,u_smc,1)));

% Minimum and maximum control 
u_excl(5) = (Lem/T0^2)*min(sqrt(dot(u_smc,u_smc,1)))*1e3;
u_excl(6) = (Lem/T0^2)*max(sqrt(dot(u_smc,u_smc,1)))*1e3;

dV_cl(2) = effort(1,3)/effort(1,1)-1;

%% Results and figures of the transfer
% Orbit representation
figure 
view(3)
hold on
plot3(Sr(1,:), Sr(2,:), Sr(3,:), 'k', 'LineWidth', 1); 
plot3(Sc_mpc(:,7), Sc_mpc(:,8), Sc_mpc(:,9), 'b', 'LineWidth', 1); 
plot3(Sc_smc(:,7), Sc_smc(:,8), Sc_smc(:,9), 'r', 'LineWidth', 1); 
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
legend('Ref.', 'MPC', 'SMC')
grid on; 
% yticklabels(strrep(yticklabels, '-', '$-$'));

figure
plot(tau, log(100*(cum_effort)/effort(1,1))); 
grid on; 
xlabel('$t$');
legend('MPC', 'SMC');
ylabel('$\Delta V(\epsilon) [\%]$');

figure_orbits = figure;
view(3)
hold on
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3),'Color','r','LineWidth', 0.9);                       % Target's orbit
plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3),'Color','b','LineWidth', 0.9);                       % Charser's initial orbit
N = plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth',1);                                                                                                  % Trasfer orbit

Q = plot3(Sc_mpc(:,1)+Sc_mpc(:,7),Sc_mpc(:,2)+Sc_mpc(:,8),Sc_mpc(:,3)+Sc_mpc(:,9),'g','LineWidth',1,'Marker','*','MarkerIndices',1:100:size(Sc_mpc,1));      % Trasfer orbit

M = plot3(Sc_smc(:,1)+Sc_smc(:,7),Sc_smc(:,2)+Sc_smc(:,8),Sc_smc(:,3)+Sc_smc(:,9),'m','LineWidth',1);                                               % Trasfer orbit

legend('Target orbit', 'Initial orbit', 'Open-loop transfer arc', 'MPC transfer arc', 'SMC transfer arc', 'Location', 'northwest', 'Autoupdate', 'off');

plot3(C(1,1),C(2,1),C(3,1),'*k');                                                                                                                   % Initial conditions
plot3(C(1,end),C(2,end),C(3,end),'*k');                                                                                                             % Final conditions
plot3(L(1,Ln), L(2,Ln), 0, '*k');
labels = {'$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$'};
text(L(1,Ln)-1e-3, L(2,Ln)-1e-3, 4e-3, labels{Ln});
hold off
grid on;

% Propulsive acceleration plot
figure;
hold on
plot(tau, sqrt(u_ol(1,:).^2+u_ol(2,:).^2+u_ol(3,:).^2)*Lem/T0^2*1e3)
plot(tau, sqrt(u_mpc(1,:).^2+u_mpc(2,:).^2+u_mpc(3,:).^2)*Lem/T0^2*1e3)
plot(tau, sqrt(u_smc(1,:).^2+u_smc(2,:).^2+u_smc(3,:).^2)*Lem/T0^2*1e3)
xlabel('$t$')
ylabel('$\|\mathbf{u}\|$')
legend('Ref.', 'MPC', 'SMC')
grid on;
% yticklabels(strrep(yticklabels, '-', '$-$'));