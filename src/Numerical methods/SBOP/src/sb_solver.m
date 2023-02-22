%% Project: 
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer using a polynomial
% shape-based approach

% Inputs: - class Problem, defining the problem of interest

% Outputs: - array C, the final state evolution matrix
%          - scalar cost, the final cost of the optimization
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar t0, the initial time of flight 
%          - scalar tf, the final time of flight 
%          - vector t, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, cost, u, t0, tf, t, exitflag, output] = sb_solver(Problem)
    % Checks 
    Problem = Problem.Check(); 

    % Setup of the algorithm
    n = Problem.PolOrder;                     % Order in the approximation of the state vector
    basis = Problem.Basis;                    % Polynomial basis to be used 
    sampling_distribution = Problem.Grid;     % Sampling grid to be used
    m = Problem.NumNodes;                     % Number of nodes in the grid
    L = Problem.DerDeg;                       % Highest derivative in the dynamics
 
    % Initial guess for the boundary control points
    mapp = 300;   
    [tapp, ~, ~, ~] = quadrature(n, mapp, sampling_distribution);

    if (~isempty(Problem.InitialGuess))
        [betaapp, t0app, tfapp, ~, Capp] = initial_approximation(Problem, basis, tapp); 
    else
        t0app = 0; 
        tfapp = 1; 
        betapp = zeros(1,1000);
    end
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(Problem, basis, tapp, Capp);
    
    % Quadrature definition
    [tau, W, ~, domain_mapping] = quadrature(n, m, sampling_distribution);

    % Final state basis
    [B, tau] = state_basis(n, L, basis, tau);

    % Initial guess reshaping
    x0 = reshape(P0, size(P0,1) * size(P0,2), []);
    x0 = [x0; t0app; tfapp; betaapp];
        
    % Objective function
    objective = @(x)cost_function(Problem, B, basis, domain_mapping, tau, W, x);

    % Non-linear constraints
    nonlcon = @(x)constraints(Problem, B, basis, domain_mapping, tau, x);

    % Upper and lower bounds 
    [P_lb, P_ub] = opt_bounds(Problem, n, size(betaapp,1));

    % Linear constraints
    [A, b, Aeq, beq] = Problem.LinConstraints(betaapp, P0);

    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, cost, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    StateCard = (max(n)+1) * Problem.StateDim;                               % Cardinal of the state modes
    P = reshape(sol(1:StateCard), Problem.StateDim, []);                     % Control points
    t0 = sol(StateCard+1);                                                   % Initial independent variable value
    tf = sol(StateCard+2);                                                   % Final independent variable value
    beta = sol(StateCard+3:end);                                             % Extra optimization parameters
    
    % Final control points imposing boundary conditions
    P = boundary_conditions(Problem, beta, t0, tf, B, basis, n, P);
    
    % Final state evolution
    C = evaluate_state(P, B, n, L);

    t = feval(domain_mapping, t0, tf, tau);                             % True domain

    % Normalization with respect to the independent variable
    m = Problem.StateDim;
    for i = 1:L
        C(1+m*i:m*(i+1),:) = C(1+m*i:m*(i+1),:) ./ (tf-t0).^i;     
    end

    u = Problem.ControlFunction(Problem.Params, beta, t0, tf, t, C);    % Control function
    t = t(1,:);
    
    % Results 
    display_results(exitflag, cost, output);
end