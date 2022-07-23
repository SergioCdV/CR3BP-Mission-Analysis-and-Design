%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer in the CR3BP using a polynomial
% shape-based approach

% Inputs: - structure system, containing the physical information of the
%           CR3BP of interest
%         - vector initial_state, the initial Cartesian state 
%         - vector final_state, the final Cartesian state
%         - scalar K, an initial desired revolutions value 
%         - scalar T, the maximum allowed acceleration
%         - scalar m, the number of sampling nodes to use 
%         - string sampling_distribution, to select the sampling distribution
%           to use 
%         - string basis, the polynomial basis to be used in the
%           optimization
%         - scalar n, the polynomial degree to be used 
%         - structure setup, containing the setup of the figures

% Outputs: - array C, the final state evolution matrix
%          - scalar dV, the final dV cost of the transfer 
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar tf, the final time of flight 
%          - scalar tfapp, the initial estimated time of flight 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_state, St, K, T, m, sampling_distribution, basis, n, setup)
    % Characteristics of the system 
    mu = system.mu;             % Characteristic gravitational parameter of the CR3BP
    t0 = system.Time;           % Fundamental time unit of the system 
    r0 = system.Distance;       % Fundamental distance unit of the system

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
    [~, Capp, Napp, tfapp] = initial_approximation(sampling_distribution, tapp, tfapp, initial, final, basis); 
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tapp, Capp, basis);
    
    % Initial guess reshaping
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    L = length(x0);
    x0 = [x0; tfapp; Napp];
    
    % Upper and lower bounds 
    P_lb = [-Inf*ones(L,1); 0; 0];
    P_ub = [Inf*ones(L,1); Inf; Inf];
    
    % Objective function
    objective = @(x)cost_function(mu, St, initial, final, n, tau, x, B, basis, sampling_distribution);
    
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Non-linear constraints
    nonlcon = @(x)constraints(mu, St, T, initial, final, n, x, B, basis, tau, sampling_distribution);
    
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
            tspan = tf*tau;
            St.Trajectory = [target_trajectory(tspan, St.Period, St.Cp); target_trajectory(tspan, St.Period, St.Cv)];
    end

    % Control input
    u = acceleration_control(mu, St, C, tf, sampling_distribution);
    u = u/tf^2;
    
    % Solution normalization
    switch (sampling_distribution)
        case 'Regularized'
            % Initial TOF 
            rapp = sqrt(Capp(1,:).^2+Capp(3,:).^2);
            tfapp = tfapp*trapz(tapp, rapp);

            % Normalised time grid
            options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
            [~, tau] = ode45(@(t,s)Sundman_transformation(basis, n, P, t, s), tau, 0, options);
    
            % Control input
            r = sqrt(C(1,:).^2+C(3,:).^2);
            u = u./(r.^2);
    
            % Final TOF 
            tf = tau(end)*tf;

        otherwise    
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
    end

    % Back transaltion of the origin of the synodic frame 
    C = cylindrical2cartesian(C, true);
    C(1:6,:) = C(1:6,:)+St.Trajectory(1:6,:);

    % Results 
    if (setup.resultsFlag)
        display_results(exitflag, output, r0, t0, tfapp, tf, dV);
    end
end
 

%% Auxiliary functions 
% Compute the derivative of time with respect to the generalized anomaly 
function [dt] = Sundman_transformation(basis, n, P, t, s)
    B = state_basis(n,s,basis);
    C = evaluate_state(P,B,n);
    dt = sqrt(C(1,:).^2+C(3,:).^2);
end
