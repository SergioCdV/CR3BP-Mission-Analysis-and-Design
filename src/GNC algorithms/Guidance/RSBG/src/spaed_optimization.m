%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer in the CR3BP using a polynomial
% shape-based approach

% Inputs: - structure system, containing the physical information of the
%           CR3BP of interest
%         - vector initial_state, the initial Cartesian state 
%         - structure St, defining the target evolution in time and the
%           vectorfield to be used in the optimization (absolute or relative
%           dynamics)
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

function [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_state, St, K, T, setup)
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

    % Target's orbit high-order approximation
    switch (St.Field)
        case 'Absolute'
            % Target state evolution
            L = [libration_points(mu) [-mu 1-mu; 0 0; 0 0; 1 0]];   % Location of the libration / primaries points
            St.Trajectory = [L(1:3,St.Center); zeros(3,1)];
            St.Trajectory = repmat(St.Trajectory,1,m);

            % Boundary conditions 
            initial = initial_state.'-St.Trajectory(:,1);           % Initial conditions
            final = St.Final.'-St.Trajectory(:,1);                  % Final conditions
            initial = cylindrical2cartesian(initial, false).';      % Initial state vector in cylindrical coordinates                  
            final = cylindrical2cartesian(final, false).';          % Final state vector in cylindrical coordinates 

        case 'Relative'
            % Target state evolution
            order = 100;                                            % Order of the polynomial

            [Cp, Cv, ~] = CTR_guidance(order, St.Trajectory(:,1).', St.Trajectory(:,2:7));

            St.Cp = Cp;                                             % Target's orbit position coordinates
            St.Cv = Cv;                                             % Target's orbit velocity coordinates
            St.Period = St.Trajectory(end,1);                       % Final target's period

            % Boundary conditions   
            initial = initial_state.'-St.Trajectory(1,2:7).';       % Relative initial conditions
            initial = cylindrical2cartesian(initial, false).';      % Initial chaser state vector in cylindrical coordinates   
            final = zeros(1,6);                                     % Rendezvous Condition

        otherwise 
            error('No valid vectorfield was chosen')
    end

    % Normlized spacecraft propulsion parameters 
    T = T*(t0^2/r0);                                        
        
    % Add additional revolutions 
    final(2) = final(2)+2*pi*K;

    % Initial guess for the boundary control points
    mapp = 300;   
    tapp = sampling_grid(mapp, sampling_distribution, '');
    [~, Capp, Napp, tfapp] = initial_approximation(mu, St, tfapp, tapp, initial, final, basis, dynamics); 
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tapp, Capp, basis);
    
    % Initial guess reshaping
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    L = length(x0);
    x0 = [x0; tfapp; Napp];
    
    % Upper and lower bounds 
    P_lb = [-Inf*ones(L,1); 0; 0];
    P_ub = [Inf*ones(L,1); 10*12*2*pi; Inf];
    
    % Objective function
    objective = @(x)cost_function(cost, mu, St, initial, final, n, tau, x, B, basis, dynamics);
    
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Non-linear constraints
    nonlcon = @(x)constraints(cost, mu, St, T, initial, final, n, x, B, basis, tau, dynamics);
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    P = reshape(sol(1:end-2), [size(P0,1) size(P0,2)]);     % Optimal control points
    tf = sol(end-1);                                        % Optimal time of flight
    N = floor(sol(end));                                    % Optimal number of revolutions 
    
    % Final control points imposing boundary conditions
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);
    
    % Final state evolution
    C = evaluate_state(P,B,n);

    % Final target trajectory 
    switch (St.Field)
        case 'Relative'
            St.Trajectory = [target_trajectory(tf, tau, St.Period, St.Cp); target_trajectory(tf, tau, St.Period, St.Cv)];
    end
    
    % Integrate the STM 
    if (setup.STM)          
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);                                            % Integration setup
        Phi0 = reshape(eye(6), [36 1]);                                                                   % Initial conditions
        [~, STM] = ode45(@(t,s)var_equations(mu, tf, St, n, P, basis, t, s), tau, Phi0, options);         % Integration
    else
        STM = [];
    end

    % Control input
    u = acceleration_control(mu, St, C, tf, dynamics);
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
            r = sundman_radius(mu, tf, St, Capp);
            tfapp = tfapp*trapz(tau, r);
    
            % Control input
            r = sundman_radius(mu, tf, St, C);
            u = u./(r.^2);
    
            % Final TOF 
            tf = tf*trapz(tau, r);

        otherwise    
    end

    % Back transaltion of the origin of the synodic frame 
    C = cylindrical2cartesian(C, true);
    C(1:6,:) = C(1:6,:)+St.Trajectory(1:6,:);
    C = [C; STM.'];

    % Results 
    if (setup.resultsFlag)
        display_results(exitflag, cost, output, r0, t0, tfapp, tf, dV);
    end
end