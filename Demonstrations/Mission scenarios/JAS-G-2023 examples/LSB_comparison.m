%% Autonomous RVD and docking in the CR3BP  %%
% Date: 06/01/23

%% Set up
set_graphics(); 
clear;
close all

%% Trajectory generation 
% Halo characteristics 
Az = 5e6:1000e3:20e6;                                                       % Orbit amplitude out of the synodic plane  
Ln = 1; 
northern_flag = 1; 

method = ["Prescribed shape-based", "Dynamics shape-based", "Dynamics shape-based", "Minimum energy", "Minimum time"];
method = method(1);

% Main computation 
[dV([1 3],:)] = cl_analysis(Az, Ln, northern_flag, method);

% Main computation 
northern_flag = -1; 
[dV([2 4],:), eps(2,:)] = performance_analysis(Az, Ln, northern_flag, method);

%%
Lem = 384400e3;                     % Mean distance from the Earth to the Moon
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);              % Normalize distances for the E-M system

%% Results and plots %% 
figure 
hold on
plot(Az, dV(1:2,:)*100);
legend('$L_1^+$','$L_1^-$')
hold off
grid on;
xlabel('$A_z$')
ylabel('$\epsilon_1 [\%]$')

figure 
hold on
plot(Az, dV(3:4,:)*100);
legend('$L_1^+$','$L_1^-$')
hold off
grid on;
xlabel('$A_z$')
ylabel('$\epsilon_2 [\%]$')

figure 
hold on
plot(Az, eps(1:2,:));
legend('$L_1^+$','$L_1^-$')
hold off
grid on; 
xlabel('$A_z$')
ylabel('$\epsilon_3$')

%% Save results
save Demonstrations\'Mission scenarios'\'JAS-G-2023 examples'\LSB_comparison.mat

%% Plots
% % Orbit representation
% figure 
% plot3(Sc(:,7), Sc(:,8), Sc(:,9), 'k', 'LineWidth', 1); 
% xlabel('$x$');
% ylabel('$y$');
% zlabel('$z$');
% grid on; 
% % yticklabels(strrep(yticklabels, '-', '$-$'));
% 
% figure_orbits = figure;
% view(3)
% hold on
% xlabel('$x$');
% ylabel('$y$');
% zlabel('$z$');
% plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3),'Color','r','LineWidth', 0.9);                    % Target's orbit
% plot3(chaser_orbit.Trajectory(:,1), chaser_orbit.Trajectory(:,2), chaser_orbit.Trajectory(:,3),'Color','b','LineWidth', 0.9);                    % Charser's initial orbit
% N = plot3(C(1,:),C(2,:),C(3,:),'k','LineWidth',1);                                                                                               % Trasfer orbit (close-loop)
% scatter3(COL(1,1:100:end),COL(2,1:100:end),COL(3,1:100:end),'k','filled','o');                                                                   % Trasfer orbit (open-loop)
% legend('Target orbit', 'Initial orbit', 'Transfer arc', 'Location', 'northwest', 'Autoupdate', 'off');
% plot3(C(1,1),C(2,1),C(3,1),'*k');                                                                                                                % Initial conditions
% plot3(C(1,end),C(2,end),C(3,end),'*k');                                                                                                          % Final conditions
% plot3(L(1,Ln), L(2,Ln), 0, '+k');
% labels = {'$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$'};
% text(L(1,Ln)-1e-3, L(2,Ln)-1e-3, 1e-2, labels{Ln});
% hold off
% grid on; 
% % yticklabels(strrep(yticklabels, '-', '$-$'));
% 
% % Propulsive acceleration plot
% figure;
% hold on
% plot(tau, sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)*Lem/T0^2*1e3)
% plot(tau, sqrt(u(4,:).^2+u(5,:).^2+u(6,:).^2)*Lem/T0^2*1e3)
% legend('OL', 'CL-SMC')
% xlabel('$t$')
% ylabel('$\|\mathbf{u}\|$')
% grid on;
% % yticklabels(strrep(yticklabels, '-', '$-$'));
% 
% figure 
% hold on
% plot(tau, rad2deg(atan2(u(2,:),u(1,:)))); 
% hold off 
% grid on;
% xlabel('Time')
% ylabel('$\theta$')
% title('Thrust in-plane angle')
% 
% figure 
% hold on
% plot(tau, rad2deg(atan2(u(3,:),sqrt(u(1,:).^2+u(2,:).^2)))); 
% hold off 
% grid on;
% xlabel('Time')
% ylabel('$\phi$')
% title('Thrust out-of-plane angle')

%% Auxiliary functions 
% Controller-assited close-loop performance
function [dV] = cl_analysis(Az, Ln, northern_flag, method_gnc)
    % Preallocation 
    dV = zeros(1,length(Az));

    % CR3BP constants 
    mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
    L = libration_points(mu);           % System libration points
    Lem = 384400e3;                     % Mean distance from the Earth to the Moon
    T0 = 28*86400;                      % Characteristic time of the Earth-Moon system
    Vc = 1.025e3;                       % Characteristic velocity of the Earth-Moon system
    n = 6;                              % Phase-space dimension
    
    % Differential corrector set up
    nodes = 10;                         % Number of nodes for the multiple shooting corrector
    maxIter = 20;                       % Maximum number of iterations
    tol = 1e-10;                        % Differential corrector tolerance

    % Main computation
    for i = 1:length(Az)
        az = dimensionalizer(Lem, 1, 1, Az(i), 'Position', 0);              % Normalize distances for the E-M system
        gamma = L(end,Ln);                                                  % Li distance to the second primary
        m = 1;                                                              % Number of periods to compute
        
        % Compute a halo seed 
        halo_param = [northern_flag az Ln gamma m];                         % Halo parameters
        [halo_seed, ~] = object_seed(mu, halo_param, 'Halo');               % Generate a halo orbit seed
        
        % Correct the seed and obtain initial conditions for a halo orbit
        [target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
        
        % Continuate the first halo orbit to locate the chaser spacecraft
        Bif_tol = 1e-2;                                                     % Bifucartion tolerance on the stability index
        num = 5;                                                            % Number of orbits to continuate
        method = 'SPC';                                                     % Type of continuation method (Single-Parameter Continuation)
        algorithm = {'Energy', NaN};                                        % Type of SPC algorithm (on period or on energy)
        object = {'Orbit', halo_seed, target_orbit.Period};                 % Object and characteristics to continuate
        corrector = 'Plane Symmetric';                                      % Differential corrector method
        direction = 1;                                                      % Direction to continuate (to the Earth)
        setup = [mu maxIter tol direction];                                 % General setup
        
        [chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
        [chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(end,:), maxIter, tol);
        
        %% Prescribed guidance using polynomial families
        % Initial conditions 
        initial_state = chaser_orbit.Trajectory(1,1:6);                     % Initial chaser conditions 
        target_state = target_orbit.Trajectory(100,1:6);                    % Initial target conditions
        
        % Minimum velocity change needed 
        Jc(1) = jacobi_constant(mu, initial_state.');
        Jc(2) = jacobi_constant(mu, target_state.');
         
        % Setup of the solution 
        T = 0.5e-3*(T0^2/Lem);                  % Maximum thrust in non-dimensional units
        GNC_LSB.Algorithm = 'Backstepping';         % Solver algorithm
        GNC_LSB.LQR.StateMatrix = 10*eye(2);        % State error weight matrix
        GNC_LSB.LQR.ControlMatrix = eye(1);         % Control effort weight matrix
        GNC_LSB.BSK.K1 = eye(2);                    % First regulation matrix
        GNC_LSB.BSK.K2 = eye(1);                    % Second regulation matrix
        GNC_LSB.Tmax = T;                           % Constrained acceleration
        GNC_LSB.TOF = pi;                           % Maneuver time
        GNC_LSB.Polynomial = 'Chebyshev';           % Polynomial family to be used
                
        % Relative solution  
        [Sr, u_ol, tf, ~] = LSB_guidance(mu, Ln, gamma, initial_state-target_state, method_gnc, GNC_LSB.TOF, GNC_LSB);  
        
        % Absolute close-loop chaser trajectory
        tau = linspace(0,tf,size(Sr,2));
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
        [~, Sc] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tau, target_state, options);
        COL = Sc.'+Sr;
        
        % Total transfer metrics 
        effort(:,1) = control_effort(tau, u_ol, false)*Vc;
                
        % Minimum and maximum control 
        u_ex(1) = (Lem/T0^2)*min(sqrt(dot(u_ol,u_ol,1)))*1e3;
        u_ex(3) = (Lem/T0^2)*max(sqrt(dot(u_ol,u_ol,1)))*1e3;
                
        % Definition of the CL controller
        GNC.Guidance.Dimension = 9;                     % Dimension of the guidance law
        GNC.Control.Dimension = 3;                      % Dimension of the control law
        GNC.Navigation.NoiseVariance = 0;               % Noise variance
        GNC.Algorithms.Guidance = '';                   % Guidance algorithm
        GNC.Algorithms.Navigation = '';                 % Navigation algorithm

        GNC.System.mu = mu;                             % Systems's reduced gravitational parameter
        GNC.System.Libration = [Ln gamma];              % Libration point ID

        GNC.Tmax = T;                                   % Maximum available acceleration
        GNC.Algorithms.Control = 'LSB';                 % Control algorithm
        GNC.LSB.Parameters = GNC_LSB;                   % Definition of the control law
        GNC.LSB.Method = method_gnc; 
        
        % Absolute close-loop chaser trajectory
        [~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tau, [target_state initial_state-target_state]);
        [~, ~, u_cl] = GNC_handler(GNC, Sc(:,1:n), Sc(:,7:end), tau); 
        
        % Total transfer metrics 
        effort(:,2) = control_effort(tau, u_cl, false)*Vc;
        
        % Minimum and maximum control 
        u_ex(2) = (Lem/T0^2)*min(sqrt(dot(u_cl,u_cl,1)))*1e3;
        u_ex(4) = (Lem/T0^2)*max(sqrt(dot(u_cl,u_cl,1)))*1e3;
                
        % Cost performance comparison
        dV(1,i) = effort(1,2)/effort(1,1)-1;
    end
end

% Controller-assited close-loop performance
function [dV, eps] = performance_analysis(Az, Ln, northern_flag, method_gnc)
    % Preallocation 
    dV = zeros(2,length(Az));
    eps = zeros(1,length(Az));

    % CR3BP constants 
    mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
    L = libration_points(mu);           % System libration points
    Lem = 384400e3;                     % Mean distance from the Earth to the Moon
    T0 = 28*86400;                      % Characteristic time of the Earth-Moon system
    Vc = 1.025e3;                       % Characteristic velocity of the Earth-Moon system
    n = 6;                              % Phase-space dimension
    
    % Differential corrector set up
    nodes = 10;                         % Number of nodes for the multiple shooting corrector
    maxIter = 20;                       % Maximum number of iterations
    tol = 1e-10;                        % Differential corrector tolerance

    % Main computation
    for i = 1:length(Az)
        az = dimensionalizer(Lem, 1, 1, Az(i), 'Position', 0);              % Normalize distances for the E-M system
        gamma = L(end,Ln);                                                  % Li distance to the second primary
        m = 1;                                                              % Number of periods to compute
        
        % Compute a halo seed 
        halo_param = [northern_flag az Ln gamma m];                         % Halo parameters
        [halo_seed, ~] = object_seed(mu, halo_param, 'Halo');               % Generate a halo orbit seed
        
        % Correct the seed and obtain initial conditions for a halo orbit
        [target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
        
        % Continuate the first halo orbit to locate the chaser spacecraft
        Bif_tol = 1e-2;                                                     % Bifucartion tolerance on the stability index
        num = 5;                                                            % Number of orbits to continuate
        method = 'SPC';                                                     % Type of continuation method (Single-Parameter Continuation)
        algorithm = {'Energy', NaN};                                        % Type of SPC algorithm (on period or on energy)
        object = {'Orbit', halo_seed, target_orbit.Period};                 % Object and characteristics to continuate
        corrector = 'Plane Symmetric';                                      % Differential corrector method
        direction = 1;                                                      % Direction to continuate (to the Earth)
        setup = [mu maxIter tol direction];                                 % General setup
        
        [chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
        [chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(end,:), maxIter, tol);
        
        %% Prescribed guidance using polynomial families
        % Initial conditions 
        initial_state = chaser_orbit.Trajectory(1,1:6);                     % Initial chaser conditions 
        target_state = target_orbit.Trajectory(100,1:6);                    % Initial target conditions
        
        % Minimum velocity change needed 
        Jc(1) = jacobi_constant(mu, initial_state.');
        Jc(2) = jacobi_constant(mu, target_state.');
        
        dVmin = velocity_norm(mu, initial_state(1:3).', Jc(1))-velocity_norm(mu, target_state(1:3).', Jc(2));
        dVmin = Vc*dVmin;
         
        % Setup of the solution 
        T = 0.5e-3*(T0^2/Lem);                  % Maximum thrust in non-dimensional units
        GNC.Algorithm = 'Backstepping';         % Solver algorithm
        GNC.LQR.StateMatrix = 10*eye(2);        % State error weight matrix
        GNC.LQR.ControlMatrix = eye(1);         % Control effort weight matrix
        GNC.BSK.K1 = eye(2);                    % First regulation matrix
        GNC.BSK.K2 = eye(1);                    % Second regulation matrix
        GNC.Tmax = T;                           % Constrained acceleration
        GNC.TOF = pi;                           % Maneuver time
        GNC.Polynomial = 'Chebyshev';           % Polynomial family to be used
                
        % Relative solution  
        [Sr, u, tf, ~] = LSB_guidance(mu, Ln, gamma, initial_state-target_state, method_gnc, GNC.TOF, GNC);  
        
        % Absolute close-loop chaser trajectory
        tau = linspace(0,tf,size(Sr,2));
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
        [~, Sc] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tau, target_state, options);
        COL = Sc.'+Sr;
        
        % Total transfer metrics 
        effort(:,1) = control_effort(tau, u, false)*Vc;
        
        % Error in time 
        [~, merit(:,1)] = figures_merit(tau, Sr);
        
        % Minimum and maximum control 
        u_ex(1) = (Lem/T0^2)*min(sqrt(dot(u,u,1)))*1e3;
        u_ex(3) = (Lem/T0^2)*max(sqrt(dot(u,u,1)))*1e3;
        
        % Regression of the guidance law 
        clear GNC; 
        order = 20; 
        [Cp, Cv, Cg] = CTR_guidance(order, tau, Sr.');
        
        GNC.Algorithms.Guidance = 'CTR';                             % Guidance algorithm
        GNC.Guidance.CTR.Order = order;                              % Order of the approximation
        GNC.Guidance.CTR.TOF = tau(end);                             % Time of flight
        GNC.Guidance.CTR.PositionCoefficients = Cp;     	         % Coefficients of the Chebyshev approximation
        GNC.Guidance.CTR.VelocityCoefficients = Cv;                  % Coefficients of the Chebyshev approximation
        GNC.Guidance.CTR.AccelerationCoefficients = Cg;              % Coefficients of the Chebyshev approximation
        GNC.Guidance.CTR.IntegralCoefficients = zeros(size(Cg));     % Coefficients of the Chebyshev approximation
        
        % Definition of the SMC controller
        GNC.Guidance.Dimension = 9;                     % Dimension of the guidance law
        GNC.Control.Dimension = 3;                      % Dimension of the control law
        GNC.Navigation.NoiseVariance = 0;               % Noise variance
        GNC.Algorithms.Guidance = 'CTR';                % Guidance algorithm
        GNC.Algorithms.Navigation = '';                 % Navigation algorithm
        GNC.System.mu = mu;                             % Systems's reduced gravitational parameter
        GNC.System.Libration = [Ln gamma];              % Libration point ID
        GNC.Tmax = T;                                   % Maximum available acceleration
        GNC.Algorithms.Control = 'SMC';                 % Control algorithm
        
        GNC.Control.SMC.Parameters = [1.000000000000000 0.985999332287318 0.006010671478548 0.1*0.013227007322678]; 
        
        % Absolute close-loop chaser trajectory
        [~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tau, [target_state initial_state-target_state], options);
        [~, ~, u(4:6,:)] = GNC_handler(GNC, Sc(:,1:n), Sc(:,7:end), tau); 
        
        % Total transfer metrics 
        effort(:,2) = control_effort(tau, u(4:6,:), false)*Vc;
        
        % Error in time 
        [e, merit(:,2)] = figures_merit(tau, Sc(:,7:12));
        
        % Minimum and maximum control 
        u_ex(2) = (Lem/T0^2)*min(sqrt(dot(u(4:6,:),u(4:6,:),1)))*1e3;
        u_ex(4) = (Lem/T0^2)*max(sqrt(dot(u(4:6,:),u(4:6,:),1)))*1e3;
        
        % Complete arc
        C = Sc(:,1:6)+Sc(:,7:12);
        C = C.';
        
        % Cost performance comparison
        dV(1,i) = effort(1,2)/effort(1,1)-1;
        dV(2,i) = effort(1,2)/dVmin-1;
    
        % Error metric comparison
        eps(1,i) = trapz(tau,sqrt(dot(Sr(1:3,:)-Sc(:,7:9).',Sr(1:3,:)-Sc(:,7:9).',1)));
    end
end