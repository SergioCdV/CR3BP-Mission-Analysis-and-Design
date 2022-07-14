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

function [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_state, final_state, K, T, m, sampling_distribution, basis, n, setup)
    % Characteristics of the system 
    mu = system.mu;             % Characteristic gravitational parameter of the CR3BP
    t0 = system.Time;           % Fundamental time unit of the system 
    r0 = system.Distance;       % Fundamental distance unit of the system

    % Approximation order 
    if (length(n) == 1)
        n = repmat(n, [1 3]);
    end

    % Boundary conditions 
    initial = cylindrical2cartesian(initial_state.', false).';      % Initial state vector in cylindrical coordinates                  
    final = cylindrical2cartesian(final_state.', false).';          % Final state vector in cylindrical coordinates 

    % Normalization
    T = T*(t0^2/r0);                                                % Spacecraft propulsion parameters 
    
    % Initial TOF
    tfapp = initial_tof(mu, T, initial_state.', final_state.');

    % Add additional revolutions 
    final(2) = final(2)+2*pi*K;

    % Initial guess for the boundary control points
    mapp = 300;   
    tapp = sampling_grid(mapp, sampling_distribution, '');
    [~, Capp, Napp, tfapp] = initial_approximation(sampling_distribution, tapp, tfapp, initial, final, basis); 
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tapp, Capp, basis);
    
    % Final collocation grid and basis 
    tau = sampling_grid(m, sampling_distribution, 'Intersection');
    [B, tau] = state_basis(n, tau, basis);

    % Initial guess reshaping
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    L = length(x0);
    x0 = [x0; tfapp; Napp];
    
    % Upper and lower bounds 
    P_lb = [-Inf*ones(L,1); 0; 0];
    P_ub = [Inf*ones(L,1); Inf; Inf];
    
    % Objective function
    objective = @(x)cost_function(mu, initial, final, n, tau, x, B, basis, sampling_distribution);
    
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Non-linear constraints
    nonlcon = @(x)constraints(mu, T, initial, final, n, x, B, basis, sampling_distribution);
    
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

    % Control input
    u = acceleration_control(mu, C, tf, sampling_distribution);
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
