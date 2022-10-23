%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 23/10/22 % 

%% Gateway Mission for MDPI %% 
% This script provides an interface to test the general control scheme for
% a rendezvous, docking and undocking mission in a re-supply mission for
% the Gateway 

% Units are non-dimensional and solutions are expressed in the 
% synodic reference frame as defined by Howell, 1984.

%% Set up %%
clc 
clear; 
close all; 

% Set up graphics 
set_graphics();

if (false)
    %% Contants and initial data %% 
    % Integration tolerances (ode113)
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    
    % Phase space dimension 
    n = 6; 
    
    % CR3BP constants 
    mu = 0.0121505;                     % Earth-Moon reduced gravitational parameter
    Lem = 384400e3;                     % Mean distance from the Earth to the Moon
    Tc = 86400*28;                      % Characteristic period of the Earth-Moon system
    Vc = 1.025e3;                       % Characteristic velocity of the Earth-Moon system
    
    % General setup 
    dt = 1e-3;                          % Time step to integrate converged trajectories
    maxIter = 20;                       % Maximum allowed iterations in the differential correction schemes
    tol = 1e-10;                        % Differential correction tolerance 
    
    %% Target and loitering halo orbit computation %%
    % Generate the L1 loitering halo 
    L = libration_points(mu);                                           % System libration points
    Az = 20e6;                                                          % Orbit amplitude out of the synodic plane. 
    Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 % Normalize distances for the E-M system
    Ln = 1;                                                             % Orbits around L1
    gamma = L(end,Ln);                                                  % Li distance to the second primary
    m = 1;                                                              % Number of periods to compute
    
    % Compute the halo orbit seed 
    halo_param = [1 Az Ln gamma m];                                     % Northern halo parameters
    [halo_seed, ~] = object_seed(mu, halo_param, 'Halo');               % Generate a halo orbit seed
    
    % Correct the seed and obtain initial conditions for a halo orbit
    [initial_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
    
    % Generate the L1 target halo
    Az = 30e6;                                                  % Orbit amplitude out of the synodic plane. Play with it! 
    Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         % Normalize distances for the E-M system
    Ln = 2;                                                     % Orbits around Li. Play with it! (L1 or L2)
    gamma = L(end,Ln);                                          % Li distance to the second primary
    m = 1;                                                      % Number of periods to compute
    param = [-1 Az Ln gamma m];                                 % Halo orbit parameters (-1 being for southern halo)
    
    % Compute the NRHO
    [halo_seed, haloT] = object_seed(mu, param, 'Halo');       
    [target_orbit, state_energy] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);
    
    %% Natural spacecraft motion %%
    r_t0 = target_orbit.Trajectory(500, 1:6);         % Initial target conditions
    r_c0 = initial_orbit.Trajectory(1, 1:6);          % Initial chaser conditions 
    rho0 = r_c0-r_t0;                                 % Initial relative conditions
    s0 = [r_t0 rho0].';                               % Initial conditions of the target and the relative state
    
    % Integration of the model
    tspan = 0:dt:2*pi; 
    [~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0, options);
    Sn = S;                
    
    % Reconstructed chaser motion 
    S_rc = S(:,1:6)+S(:,7:12);                        % Reconstructed chaser motion via Encke method
    
    % Insertion plot 
    figure
    view(3)
    hold on
    plot3(Sn(:,1), Sn(:,2), Sn(:,3));
    plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3));
    scatter3(Sn(1,1), Sn(1,2), Sn(1,3), 'blue', 'filled');
    scatter3(S_rc(1,1), S_rc(1,2), S_rc(1,3), 'red', 'filled');
    hold off
    grid on;
    legend('Target orbit', 'Loiter orbit', 'Target at HOI', 'Chaser at HOI'); 
    
    %% Initial phase: SMC + continuous iLQR 
    % Guidance computation
    GNC.Guidance.Dimension = 9;                     % Dimension of the guidance law
    GNC.Control.Dimension = 3;                      % Dimension of the control law
    GNC.Navigation.NoiseVariance = 0;               % Noise variance
    
    GNC.Algorithms.Guidance = '';                   % Guidance algorithm
    GNC.Algorithms.Navigation = '';                 % Navigation algorithm
    
    GNC.System.mu = mu;                             % Systems's reduced gravitational parameter
    GNC.System.Libration = [Ln gamma];              % Libration point ID
    
    GNC.Algorithms.Control = 'SDRE';                % Control algorithm
    GNC.Control.SDRE.Q = 5e2*blkdiag(eye(3), 1e-6*eye(6));                          % Control algorithm
    GNC.Control.SDRE.M = eye(3);                    % Penalty on the control effort
    GNC.Control.SDRE.Model = 'RLM';                 % Control algorithm model
    
    GNC.Tmax = 3e-3 / (4*pi^2*Lem/Tc^2);            % Maximum available acceleration
    
    % Augmented initial conditions 
    int = zeros(3,1);                               % Integral of the relative position
    s0 = [s0; int].';                               % Initial conditions
    
    dt = 1e-3;                                      % Zero-order hold time span
    tf(1) = pi;                                     % Approach TOF
    tspan_1 = 0:dt:tf(1);
    
    [~, S1] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan_1, s0, options); 
    
    % % Offline Guidance regression
    % order = 50; 
    % [Cp, Cv, Cg] = CTR_guidance(order, tspan_1, Sg(:,n+1:2*n));
    % 
    % GNC.Algorithms.Guidance = 'CTR';                             % Guidance algorithm
    % GNC.Guidance.CTR.Order = order;                              % Order of the approximation
    % GNC.Guidance.CTR.TOF = tspan(end);                           % Time of flight
    % GNC.Guidance.CTR.PositionCoefficients = Cp;     	         % Coefficients of the Chebyshev approximation
    % GNC.Guidance.CTR.VelocityCoefficients = Cv;                  % Coefficients of the Chebyshev approximation
    % GNC.Guidance.CTR.AccelerationCoefficients = Cg;              % Coefficients of the Chebyshev approximation
    % GNC.Guidance.CTR.IntegralCoefficients = zeros(size(Cg));     % Coefficients of the Chebyshev approximation
    % 
    % % Control 
    % GNC.Navigation.NoiseVariance = 1/Lem;
    % 
    % GNC.Algorithms.Control = 'SMC';                              % Control algorithm
    % 
    % GNC.Control.SMC.Parameters = [1.000000000000000 0.985999332287318 0.006010671478548 0.1*0.013227007322678]; 
    % 
    % tic;
    % [~, S1] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan_1, s0, options);
    % Time(1) = toc; 
    
    % Control law
    [~, ~, u_smc] = GNC_handler(GNC, S1(:,1:n), S1(:,n+1:end), tspan_1);
    
    % Error in time 
    [e_smc, merit_smc] = figures_merit(tspan_1, S1(:,n+1:2*n));
    
    % Control integrals
    effort_smc = control_effort(tspan_1, u_smc, false);
    
    % Range reduction 
    range(1) = 1-norm(S1(end,n+1:n+3))/norm(S1(1,n+1:n+3));
    
    figure
    view(3)
    hold on
    plot3(Sn(:,1), Sn(:,2), Sn(:,3));
    plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3));
    plot3(S1(:,1)+S1(:,7), S1(:,2)+S1(:,8), S1(:,3)+S1(:,9));
    legend('Target orbit', 'Loiter orbit', 'SDRE transfer', 'AutoUpdate', 'off'); 
    % plot3(Sg(:,1)+Sg(:,7), Sg(:,2)+Sg(:,8), Sg(:,3)+Sg(:,9));
    scatter3(Sn(1,1), Sn(1,2), Sn(1,3), 'blue', 'filled');
    scatter3(S_rc(1,1), S_rc(1,2), S_rc(1,3), 'red', 'filled');
    scatter3(L(1,1), L(2,1), 0, 'k', 'filled')
    scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
    hold off
    text(L(1,1)+1e-2, L(2,Ln), 1e-4, '$L_1$');
    text(L(1,Ln)-2e-2, L(2,Ln), 1e-4, '$L_2$');
    grid on;
    
    figure 
    view(3)
    hold on 
    % plotObject(Sg(:,7:12)); 
    plotObject(S1(:,7:12));
    
    %% Second phase: MPC + Opt-MI
    % Set up of the optimization
    method = 'NPL';                               % Method to solve the problem
    core = 'Linear';                              % Number of impulses
    cost_function = 'Position';                   % Target a position rendezvous
    
    tf(2) = 1;                                  % Time of flight
    
    % Thruster characteristics 
    Tmin = -1e-3;                                 % Minimum thrust capability (in velocity impulse)
    Tmax = 1e-3;                                  % Maximum thrust capability (in velocity impulse)
    
    % Main computation 
    s0 = S1(end,1:2*n);
    tic
    [tspan_21, S21, u_mpc1, state{2}] = MPC_control(mu, cost_function, Tmin, Tmax, tf(2), s0, core, method);
    Time(2) = toc;
    
    s0 = S21(end,1:2*n);
    tic
    [tspan_22, S22, u_mpc2, state{2}] = MPC_control(mu, cost_function, Tmin, Tmax, 0.25, s0, core, method);
    Time(2) = Time(2) + toc;
    
    % Complete trajectory 
    S2 = [S21(1:end-1,1:2*n); S22(1:end,1:2*n)];
    tspan_2 = [tspan_21(1:end-1) tspan_21(end)+tspan_22];
    tf(2) = tf(2) + 0.25;
    u_mpc = [u_mpc1 u_mpc2];
    
    % Error in time 
    [e_mpc, merit_mpc] = figures_merit(tspan_2, S2(:,n+1:2*n));
    
    % Control integrals
    effort_mpc = control_effort(tspan_2, u_mpc, true);
    
    % Range reduction 
    range(2) = 1-norm(S2(end,n+1:n+3))/norm(S2(1,n+1:n+3));
    
    %% Third phase: formation-flying
    % Guidance (Lissajous trajectory around the target)
    tf(3) = pi; 
    
    % Generate an 1 km Lissajous relative orbit
    tspan_3 = (0:1e-3:tf(3)).';
    Ax = 1e3/Lem;
    Az = 1e3/Lem;
    Sg(:,1:3) = [Ax*cos((2*pi/0.5)*tspan_3) Ax*sin((2*pi/0.5)*tspan_3) Az*sin((2*pi/3)*tspan_3)];
    Sg(:,4:6) = [-(2*pi/0.5)*Ax*sin((2*pi/0.5)*tspan_3) (2*pi/0.5)*Ax*cos((2*pi/0.5)*tspan_3) (2*pi/3)*Az*cos((2*pi/3)*tspan_3)];
    
    % Compute the trajectory as a Chebyshev analytical expression
    tspan_3 = tspan_3.';
    order = 50; 
    [Cp, Cv, Cg] = CTR_guidance(order, tspan_3, Sg);
    
    % GNC algorithms definition 
    GNC.Algorithms.Guidance = 'CTR';            % Guidance algorithm
    GNC.Algorithms.Navigation = '';             % Navigation algorithm
    GNC.Algorithms.Control = 'SMC';             % Control algorithm
    FNC.Algorithms.Solver = 'Encke';            % Dynamics vector field to be solved
    GNC.Guidance.Dimension = 9;                 % Dimension of the guidance law
    GNC.Control.Dimension = 3;                  % Dimension of the control law
    GNC.System.mu = mu;                         % System reduced gravitational parameter
    
    GNC.Navigation.NoiseVariance = dimensionalizer(Lem, 1, 1, 1, 'Position', 0);
    
    % Controller parameters
    GNC.Control.SMC.Parameters = [1.000000000000000 0.432562054680836 0.070603623964497 0.099843662546135]; 
    
    % Guidance parameters 
    GNC.Guidance.CTR.Order = order;                     % Order of the approximation
    GNC.Guidance.CTR.TOF = tf(3);                       % Time of flight
    GNC.Guidance.CTR.PositionCoefficients = Cp;     	% Coefficients of the Chebyshev approximation
    GNC.Guidance.CTR.VelocityCoefficients = Cv;         % Coefficients of the Chebyshev approximation
    GNC.Guidance.CTR.AccelerationCoefficients = Cg;     % Coefficients of the Chebyshev approximation
    
    % Controller
    s0 = S2(end,1:2*n); 
    tic
    [~, S3] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan_3, s0, options);
    Time(3) = toc;
    
    % Control effort
    [~, ~, u] = GNC_handler(GNC, S3(:,1:n), S3(:,n+1:2*n), tspan_3); 
    effort_ff = control_effort(tspan_3, u, false);
    
    %% Final results 
    % Complete trajectory 
    St = [S1(1:end-1,1:2*n); S2(1:end-1,1:2*n); S3(1:end,1:2*n)];
    
    % Total integration time
    tspan = [tspan_1(1:end-1) tspan_1(end)+tspan_2(1:end-1) tspan_2(end)+tspan_3(1:end)];
    
    save Gateway_mission;
else
    load Gateway_mission;
end

%% Plotting
figure 
view(3) 
hold on 
t = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'k*', 'MarkerIndices', 1:200:size(Sn,1));
c = plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3), 'k-*', 'MarkerIndices', 1:200:size(S_rc,1));
plot3(St(length(tspan_1):end,1)+St(length(tspan_1):end,7), St(length(tspan_1):end,2)+St(length(tspan_1):end,8), St(length(tspan_1):end,3)+St(length(tspan_1):end,9), 'r', 'LineWidth', 0.7); 
legend('Target orbit', 'Chaser motion', 'Location', 'northwest', 'AutoUpdate', 'off');
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
hold off
text(L(1,Ln)-1e-4, L(2,Ln), 1e-4, '$L_2$');
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;

figure 
view(3) 
plot3(St(:,7), St(:,8), St(:,9));
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
grid on;

% Configuration space evolution (Phase I)
figure
subplot(1,2,1)
hold on
plot(tspan_1, S1(:,7)); 
plot(tspan_1, S1(:,8)); 
plot(tspan_1, S1(:,9)); 
hold off
xlabel('$t$');
ylabel('$\mathbf{\rho}$');
grid on;
legend('$x$', '$y$', '$z$');
subplot(1,2,2)
hold on
plot(tspan_1, S1(:,10)); 
plot(tspan_1, S1(:,11)); 
plot(tspan_1, S1(:,12)); 
hold off
xlabel('$t$');
ylabel('$\dot{\mathbf{\rho}}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');

% Configuration space evolution (Phase II)
figure
subplot(1,2,1)
hold on
plot(tspan_2(1:120), S2(1:120,7)); 
plot(tspan_2(1:120), S2(1:120,8)); 
plot(tspan_2(1:120), S2(1:120,9)); 
hold off
xlabel('$t$');
ylabel('$\mathbf{\rho}$');
grid on;
legend('$x$', '$y$', '$z$');
subplot(1,2,2)
hold on
plot(tspan_2(1:120), S2(1:120,10)); 
plot(tspan_2(1:120), S2(1:120,11)); 
plot(tspan_2(1:120), S2(1:120,12)); 
hold off
xlabel('$t$');
ylabel('$\dot{\mathbf{\rho}}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');

% Configuration space evolution (Phase III)
figure
subplot(1,2,1)
hold on
plot(tspan_3, S3(:,7)); 
plot(tspan_3, S3(:,8)); 
plot(tspan_3, S3(:,9)); 
hold off
xlabel('$t$');
ylabel('$\mathbf{\rho}$');
grid on;
legend('$x$', '$y$', '$z$');
subplot(1,2,2)
hold on
plot(tspan_3, S3(:,10)); 
plot(tspan_3, S3(:,11)); 
plot(tspan_3, S3(:,12));
hold off
xlabel('$t$');
ylabel('$\dot{\mathbf{\rho}}$');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');

if (false)
    dh = 50; 
    steps = fix(size(St,1)/dh);
    M = cell(1,steps);
    h = figure;
    filename = 'webb.gif';
    view([40 20])
    hold on
    plot3(flip(St0(:,1)), flip(St0(:,2)), flip(St0(:,3)), '.r', 'Linewidth', 0.1);
    plot3(St(:,1), St(:,2), St(:,3), '.-b', 'Linewidth', 0.1);
    plot3(St4(:,1)+St4(:,7), St4(:,2)+St4(:,8), St4(:,3)+St4(:,9), '.r', 'Linewidth', 0.1); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
    scatter3(1-mu, 0, 0, 'k', 'filled');
    text(L(1,Ln)-2e-3, L(2,Ln), 0, '$L_2$');
    text(1-mu-1e-3, 0, 1e-3, '$M_2$');
    grid on;
    title('Rendezvous simulation');
    
    for i = 1:dh:size(St,1)
        T = scatter3(St(i,1), St(i,2), St(i,3), 20, 'b', 'filled'); 
        V = scatter3(St(i,1)+St(i,7), St(i,2)+St(i,8), St(i,3)+St(i,9), 20, 'r', 'filled');
        drawnow;
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); 
        if (i == 1) 
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1e-3); 
        else 
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1e-3); 
        end 
        delete(T); 
        delete(V);
    end
    hold off
end
