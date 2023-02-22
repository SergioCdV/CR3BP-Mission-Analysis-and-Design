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
Az = 20e6:100e4:30e6;                    % Orbit amplitude out of the synodic plane 
Ln = 2;                                 % Libration point index
northern_flag = -1;                      % Northern/southern familiy flag
method = 'Minimum energy';              % Direct prescritpion 

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

% Preallocation 
dV = zeros(1,length(Az));
u_ex = zeros(4,length(Az));
sigma = zeros(2,length(Az));
epsilon = zeros(1,length(Az));
Jc = zeros(2,length(Az)); 

% Target orbit
az = dimensionalizer(Lem, 1, 1, Az(1), 'Position', 0);              % Normalize distances for the E-M system
halo_param = [northern_flag az Ln gamma m];                         % Halo parameters
[halo_seed, ~] = object_seed(mu, halo_param, 'Halo');               % Generate the halo orbit seed

% Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

target_state = target_orbit.Trajectory(100,1:6);                    % Initial target conditions
Jc(2,:) = repmat(jacobi_constant(mu, target_state.'),1,length(Az)); % Jacobi constant

options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);

% Main computation 
for i = 1:length(Az)
    az = Az(i)/Lem;                                                     % Normalize distances for the E-M system
    
    % Compute a halo seed 
    halo_param = [northern_flag az Ln gamma m];                         % Halo parameters
    [halo_seed, ~] = object_seed(mu, halo_param, 'Halo');               % Generate a halo orbit seed
    
    % Correct the seed and obtain initial conditions for a halo orbit
    [chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
     
    initial_state = chaser_orbit.Trajectory(1,1:6);                     % Initial chaser conditions 
    
    % Minimum velocity change needed 
    Jc(1,i) = jacobi_constant(mu, initial_state.');
                     
    % Relative solution  
    rho = initial_state-target_state;
    sigma(1,i) = norm(rho(1:3));
    sigma(2,i) = norm(rho(4:6));
    [Sr, u_ol, tf, ~] = LSB_guidance(mu, Ln, gamma, rho, method, GNC_LSB.TOF, GNC_LSB);  
    
    % Absolute target trajectory
    tau = linspace(0,tf,size(Sr,2));
    [~, Sc] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tau, target_state, options);
    
    % Total transfer metrics 
    effort(:,1) = control_effort(tau, u_ol, false) * Vc;
            
    % Minimum and maximum control 
    u_ex(1,i) = (Lem/T0^2)*min(sqrt(dot(u_ol,u_ol,1)))*1e3;
    u_ex(2,i) = (Lem/T0^2)*max(sqrt(dot(u_ol,u_ol,1)))*1e3;
                    
    u_smc = zeros(3,length(tau));
    for j = 1:length(tau)
        aux = nlr_model(mu, true, false, false, 'Encke', tau, [Sc(j,:).'; Sr(1:n,j)]);
        u_smc(:,j) = Sr(7:9,j)-aux(10:12);
    end

    % Minimum and maximum control 
    u_ex(3,i) = (Lem/T0^2)*min(sqrt(dot(u_smc,u_smc,1)))*1e3;
    u_ex(4,i) = (Lem/T0^2)*max(sqrt(dot(u_smc,u_smc,1)))*1e3;

    % Feasibility analysis 
    feas = u_ol-u_smc;
    epsilon(1,i) = Vc * trapz(tau, sqrt(dot(feas(1:3,:),feas(1:3,:),1)));
    dV(1,i) = epsilon(1,i)/effort(1,1)-1;
end

%% Results and figures 
% Feasibility plot 
figure 
plot(sigma(1,:), epsilon);
grid on; 
ylabel('$\delta \Delta V$')
xlabel('$\|\rho_0\|$')

figure 
plot(sigma(1,:), dV*100);
grid on; 
ylabel('$\delta \Delta V$ [\%]')
xlabel('$\|\rho_0\|$')

figure 
plot(sigma(1,:), u_ex([1 3],:));
grid on; 
legend('$u_{A,min}$', '$u_{E,min}$')
xlabel('$\|\rho_0\|$')
ylabel('$\|\mathbf{u}\|$')

figure 
plot(sigma(1,:), u_ex([2 4],:));
grid on;
ylabel('$\|\mathbf{u}\|$')
legend('$u_{A,max}$', '$u_{E,max}$')
xlabel('$\|\rho_0\|$')