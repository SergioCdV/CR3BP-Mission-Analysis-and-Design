%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/08/22
% File: yrgb_guidance.m 
% Issue: 0 
% Validated: 20/08/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer in the CR3BP using a polynomial
% shape-based approach for halo raising usign hyperbolic manifolds

% Inputs: - structure system, containing the physical information of the
%           CR3BP of interest
%         - vector initial_state, the initial Cartesian state 
%         - structure Sc, defining the reference orbit evolution in time and the
%           vectorfield to be used in the optimization
%         - scalar K, an initial desired revolutions value 
%         - scalar T, the maximum allowed acceleration
%         - structure setup, containing the setup of the algorithm in general

% Outputs: - array C, the final state evolution matrix
%          - scalar dV, the final dV cost of the transfer 
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar tf, the final time of flight 
%          - scalar tfapp, the initial estimated time of flight 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, dV, u, tf, tfapp, tau, exitflag, output] = yrsb_optimization(system, curve, Sc, K, T, setup)
    % Setup of the algorithm
    n = setup.order;                        % Order in the approximation of the state vector
    dynamics = setup.formulation;           % Formulation of the dynamics
    basis = setup.basis;                    % Polynomial basis to be used 
    sampling_distribution = setup.grid;     % Sampling grid to be used
    m = setup.nodes;                        % Number of nodes in the grid
    cost = setup.cost_function;             % Cost function to be minimized   

    % Characteristics of the system 
    mu = system.mu;                         % Characteristic gravitational parameter of the CR3BP
    t0 = system.Time;                       % Fundamental time unit of the system 
    r0 = system.Distance;                   % Fundamental distance unit of the system

    % Approximation order 
    if (length(n) == 1)
        n = repmat(n, [1 3]);
    end

    % Initial TOF
    tfapp = 2*pi;

    % Final collocation grid and basis 
    tau = sampling_grid(m, sampling_distribution, '');
    [B, tau] = state_basis(n, tau, basis);

    % Target state evolution
    order = 20;                                            % Order of the polynomial

    theta = linspace(0,2*pi,size(Sc.Trajectory,1));
    [Cp, Cv, ~] = CTR_guidance(order, theta, Sc.Trajectory(:,2:7));

    Sc.Cp = Cp;                                             % Target's orbit position coordinates
    Sc.Cv = Cv;                                             % Target's orbit velocity coordinates
    Sc.Period = Sc.Trajectory(end,1);                       % Final target's period

    % Compute the initial unstable vector displacement 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  % Integration setup 
    Sc0 = [Sc.Trajectory(1,2:7) reshape(eye(6), [1 36])];   % Chaser orbit initial conditions
    tspan = 0:1e-3:Sc.Period;                               % Integration time span

    [~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, Sc0, options);
 
    M = reshape(S(end,7:end), [6 6]);                       % Monodromy matrix 
    [V, ~] = eig(M);                                        % Eigenspectrum of the monodromy matrix
    us = V(:,1)/norm(V(:,1));                               % Unstable manifold
    epsilon = 1e-4;                                         % Unstable manifold displacement

    if (Sc.Branch == 'L')
        epsilon = -epsilon; 
    end

    % Boundary conditions   
    initial = epsilon*us.';                                                                                % Relative initial conditions in cartesian coordinates 
    final = final_orbit(curve, 0);                                                                         % Final initial guess for the insertion point
    final = final-target_trajectory(sampling_distribution, tfapp, tau(end), Sc.Period, [Sc.Cp; Sc.Cv]);    % Final initial guess for the insertion point
    final = cylindrical2cartesian(final, false).';                                                         % Final conditions in cylindrical coordinates 

    % Normlized spacecraft propulsion parameters 
    T = T*(t0^2/r0);                                        
        
    % Add additional revolutions 
    final(2) = final(2)+2*pi*K;

    % Initial guess for the boundary control points
    mapp = 300;   
    tapp = sampling_grid(mapp, sampling_distribution, '');
    [~, Capp, Napp, tfapp] = initial_approximation(mu, Sc, tfapp, tapp, initial, final, basis, dynamics); 
     
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tapp, Capp, basis);
    
    % Initial guess reshaping
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    L = length(x0);
    x0 = [x0; tfapp; Napp; pi/2];
    
    % Upper and lower bounds 
    P_lb = [-Inf*ones(L,1); 0; 0; 0];
    P_ub = [Inf*ones(L,1); 10*12*2*pi; Inf; 2*pi];
    
    % Objective function
    objective = @(x)cost_function(curve, cost, mu, Sc, initial, n, tau, x, B, basis, sampling_distribution, dynamics);
    
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Non-linear constraints
    nonlcon = @(x)constraints(curve, cost, mu, Sc, T, initial, n, x, B, basis, tau, sampling_distribution, dynamics);
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    P = reshape(sol(1:end-3), [size(P0,1) size(P0,2)]);     % Optimal control points
    tf = sol(end-2);                                        % Optimal time of flight on the unstable manifold
    N = floor(sol(end-1));                                  % Optimal number of revolutions 
    theta = sol(end);                                       % Optimal insertion phase

    % Final target trajectory 
    Sc.Trajectory = target_trajectory(sampling_distribution, sum(tf), tau, Sc.Period, [Sc.Cp; Sc.Cv]);

    % Compute the insertion phase and final conditions
    final = final_orbit(curve, theta);                      % Final initial guess for the insertion point
    final = final-Sc.Trajectory(:,end);                     % Final initial guess for the insertion point
    final = cylindrical2cartesian(final, false).';
    
    % Final control points imposing boundary conditions
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);
    
    % Final state evolution
    C = evaluate_state(P,B,n);
    
    % Integrate the STM 
    if (setup.STM)          
        STM = stm_computation(mu, tf, Sc, n, P, B, sampling_distribution, basis, tau, 'Numerical'); 
    else
        STM = [];
    end

    % Control input
    u = acceleration_control(mu, Sc, C, tf, dynamics);
    u = u/tf^2;

    % Time domain normalization 
    switch (sampling_distribution)
        case 'Chebyshev'
            tau = (1/2)*(1+tau);
            tf = tf*2;
        case 'Legendre'
            tau = (1/2)*(1+tau);
            tf = tf*2;
        case 'Laguerre'
            tau = collocation_grid(m, 'Legendre', '');
            tau = (1/2)*(1+tau);
            tf = tf*2;
    end
    
    % Solution normalization
    switch (dynamics)
        case 'Sundman'
            % Initial TOF 
            r = sundman_radius(mu, tf, Sc, Capp);
            tfapp = tfapp*trapz(tau, r);
    
            % Control input
            r = sundman_radius(mu, tf, Sc, C);
            u = u./(r.^2);
    
            % Final TOF 
            tf = tf*trapz(tau, r);

        otherwise    
    end

    % Back transaltion of the origin of the synodic frame 
    C = cylindrical2cartesian(C, true);
    C(1:6,:) = C(1:6,:)+Sc.Trajectory(1:6,:);
    C = [C; STM];

    % Results 
    if (setup.resultsFlag)
        display_results(exitflag, cost, output, r0, t0, tfapp, sum(tf), dV);
    end
end