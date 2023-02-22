%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Cost function %%
% Function to compute the cost function to be minimized

% Inputs: - class Problem, defining the problem at hands
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - function handle domain_mapping, to map the independent variable
%           domain to the original one
%         - vector tau, the vector of collocation points
%         - vector W, the quadrature weights
%         - vector x, the degree of freedom to be optimized
%         - string dynamics, the independent variable parametrization to be
%           used

% Outputs: - scalar r, the cost index to be optimized

function [r] = cost_function(Problem, B, basis, domain_mapping, tau, W, x)
    % Optimization variables
    L = Problem.DerDeg;                                                 % Maximum derivative degree
    n = Problem.PolOrder;                                               % Order of the polynomial approximation
    m = Problem.StateDim;                                               % State dimension
    StateCard = (max(n)+1) * m;                                         % Cardinal of the state modes
    P = reshape(x(1:StateCard), m, []);                                 % Control points
    t0 = x(StateCard+1);                                                % Initial independent variable value
    tf = x(StateCard+2);                                                % Final independent variable value
    beta = x(StateCard+3:end);                                          % Extra optimization parameters
    
    % Evaluate the boundary conditions
    P = boundary_conditions(Problem, beta, t0, tf, B, basis, n, P);     % Boundary conditions control points
    s = evaluate_state(P, B, n, L);                                     % State evolution
    t = feval(domain_mapping, t0, tf, tau);                             % Original time independent variable

    % Normalization
    for i = 1:L
        s(1+m*i:m*(i+1),:) = s(1+m*i:m*(i+1),:) ./ (tf-t0)^i;     
    end

    u = Problem.ControlFunction(Problem.Params, beta, t0, tf, t, s);    % Control function
        
    % Evaluate the cost function (Lagrange and Mayer terms)
    [M, L] = Problem.CostFunction(Problem.Params, beta, t0, tf, s, u); 

    if (isempty(W))
        r = M + (tf-t0) * trapz(tau,L);
    else
        r = M + (tf-t0) * dot(W,L);
    end
end
