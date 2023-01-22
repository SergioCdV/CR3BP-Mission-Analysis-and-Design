%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 15/01/23 % 

%% Relative Floquet Stationkeeping demonstration %% 
% This script provides an interface to test the RFSK strategies for optimal long term stationkeeping

% Units are non-dimensional and solutions are expressed in the synodic
% reference frame as defined by Howell, 1984.

%% Set up %%
clc
clear; 
close all; 

% Set up graphics 
set_graphics();

% Integration tolerances (ode113)
if (0)
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
    [~, Sn] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);   
    
    % Guidance law for the Jacobi constant
    Jref = jacobi_constant(mu, Sn(1,1:n).');
    
    %% GNC-RFSK algorithm definition 
    GNC.Algorithms.Guidance = '';                         % Guidance algorithm
    GNC.Algorithms.Navigation = '';                       % Navigation algorithm
    GNC.Algorithms.Control = 'RFSK';                      % Control algorithm
    
    GNC.Guidance.Dimension = 9;                           % Dimension of the guidance law
    GNC.Control.Dimension = 3;                            % Dimension of the control law
    
    GNC.Navigation.NoiseVariance = 0;
    
    GNC.System.mu = mu;                                   % Systems's reduced gravitational parameter
    GNC.Control.RFSK.method = 'LQR';                      % Solver to be used
    GNC.Control.RFSK.Q = 1*eye(2);                        % Penalty on the state error
    GNC.Control.RFSK.M = 5e-3*eye(3);                     % Penalty on the control effort
    GNC.Control.RFSK.Reference = Jref;                      % Reference state
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
    
    %% GNC-LQR algorith definition
    GNC.System.Libration = [Ln gamma];                                 % Libration point ID
    GNC.Control.LQR.Model = 'RLM';                                     % LQR model
    GNC.Control.LQR.Q = 2*eye(9);                                      % Penalty on the state error
    GNC.Control.LQR.M = 2e-2*eye(3);                                   % Penalty on the control effort
    GNC.Control.LQR.Reference = target_orbit.Trajectory(1,1:6);        % Reference operting point
    
    %% Montecarlo analysis
    % Initial conditions 
    k = [100e3*ones(1,3)/Lem 0.5*ones(1,3)/Vc];
    
    r_t0 = Sn(1,1:6);                                         % Initial guidance target conditions
    tspan = 0:dt:tf;                                          % New integration time span
    
    % Compute the trajectory
    iter = 1e3; 
    Error = zeros(1,iter);      % Initial error
    Time = zeros(2,iter);       % Computational cost
    Merit = zeros(4,iter);      % Error cost
    Effort = zeros(2,iter);     % Control cost
    
    % Nominal trajectory
    s0 = [r_t0 zeros(1,6) reshape(eye(n), [1 n^2])];
    s0_lqr = [s0(1:12) zeros(1,3) s0(13:end)]; 
    [~, Sr] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0(1:2*n), options);
    
    K = 0.1;
    
    % Monte Carlos analysis
    for i = 1:iter
        % Initial insertion error
        s0(n+1:2*n) = normrnd(zeros(1,6),k);     % Noisy relative initial conditions
        s0_lqr(n+1:2*n) = s0(n+1:2*n);           % Noisy relative initial conditions
        Error(1,i) = norm(s0(7:9));
        Error(2,i) = norm(s0(10:12));
        
        % Stationkeeping trajectory. RFSK
        GNC.Algorithms.Control = 'RFSK';
        tic
        [~, Staux1] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan(tspan < K*target_orbit.Period), s0, options);  
        [~, Staux2] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan(tspan >= K*target_orbit.Period), Staux1(end,:), options);
        Time(1,i) = toc;
        
        % Error in time 
        Merit(1:2,i) = [norm(Staux2(end,7:9)); norm(Staux2(end,10:12));];
        
        % Control law
        [~, ~, u] = GNC_handler(GNC, Staux1(:,1:n), Staux1(:,n+1:end), tspan(tspan < K*target_orbit.Period));
        
        % Control integrals
        effort = control_effort(tspan(tspan < K*target_orbit.Period), u, false)*Vc;
        Effort(1,i) = effort(1);
    
        % Stationkeeping trajectory. LQR
        GNC.Algorithms.Control = 'LQR';
        tic
        [~, Staux1] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan(tspan < K*target_orbit.Period), s0_lqr(1:15), options);  
        [~, Staux2] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan(tspan >= K*target_orbit.Period), Staux1(end,:), options);
        Time(2,i) = toc;
        
        % Error in time 
        Merit(3:4,i) = [norm(Staux2(end,7:9)); norm(Staux2(end,10:12));];
        
        % Control law
        [~, ~, u] = GNC_handler(GNC, Staux1(:,1:n), Staux1(:,n+1:end), tspan(tspan < K*target_orbit.Period));
        
        % Control integrals
        effort = control_effort(tspan(tspan < K*target_orbit.Period), u, false)*Vc;
        Effort(2,i) = effort(1);
    end
    
    save MonteCarlo;
else
    load MonteCarlo;
end

%% Results %% 
% Ordering 
[Error, index] = sort(Error(1,:)); 
Effort = Effort(:,index);
Merit = Merit(:,index);
Time = Time(:,index);

% Control cost figure 
figure 
hold on
scatter(Error(1,:), Effort(1,:), 'filled', 'b')
scatter(Error(1,:), Effort(2,:), 'filled', 'r')
legend('RFSK', 'LQR'); 
grid on; 
xlabel('Initial velocity error $\dot{\mathbf{\rho}_0}$')
ylabel('$\Delta V$ [m/s]')

% Error cost figure 
figure 
hold on
scatter(Error(1,:), log(Merit(1,:)), 'filled')
scatter(Error(1,:), log(Merit(2,:)), 'filled')
scatter(Error, log(Merit(3,:)), '*')
scatter(Error, log(Merit(4,:)), '*')
legend('RFSK-$\rho$', 'RFSK-$\dot{\rho}$', 'LQR-$\rho$', 'LQR-$\dot{\rho}$'); 
grid on; 
xlabel('Initial velocity error $\dot{\mathbf{\rho}_0}$')
ylabel('Error log $e(t_f)$')

% Computational cost figure
figure 
hold on
scatter(Error(1,:), Time(1,:), 'filled', 'b')
scatter(Error(1,:), Time(2,:), 'filled', 'r')
legend('RFSK', 'LQR'); 
grid on; 
xlabel('Initial velocity error $\dot{\mathbf{\rho}_0}$')
ylabel('Computational cost $T_c$ [s]')
