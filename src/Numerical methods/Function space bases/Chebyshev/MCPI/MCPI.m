%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 28/09/22
% File: MCPI.m 
% Issue: 0 
% Validated: 

%% Modified Chebyshev Picard Iterations %%
% This function allows to compute a given quadrature using Bai's vectorized
% MCPI method

% Inputs: - string kind, specifying the kind of Chebyshev polynomials to
%           compute
%         - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated

% Outpus: - vector Pn, containing the evaluated Chebyshev polynomials 

order = 100; 
tau = flip(cos((0:order)*pi/order));

x0 = [1000 -3]; 

t = (1/2)*(tau+1);
[t,x1] = ode45(@(t,x)([-2*t*x(1); -x(2)]), t, x0, odeset('RelTol',2.25E-14, 'AbsTol', 1e-22));


dynamics = @(tau,x)([-x(:,1).*(1/2*(tau+1)) -x(:,2)]);
[tau, x, state] = MCPI_IVP(tau, repmat(x0,order+1,1), dynamics, order, 1e-22);

function [tau, x, state] = MCPI_IVP(tau, x0, dynamics, order, tol)
    % Set up of the method 
    GoOn = true;                            % Convergence boolean
    iter = 1;                               % Initial iteration
    maxIter = 1e3;                          % Maximum number of iterations
    error = 1e3;                            % Initial error

    % Constants 
    N = order+1;                            % Number of approximation terms

    W = diag(ones(1,N));                    % Weight matrix for the state vector
    W(1,1) = 1/2;                           % Weight matrix for the state vector

    T = CH_basis('first', order, tau).';    % Chebyshev basis

    V = (2/N)*diag(ones(1,N));              % Approximation weights for the dynamics
    V(1,1) = 1/N;                           % Approximation weights for the dynamics
    V(N,N) = 1/N;                           % Approximation weights for the dynamics

    R = (1/2)./(1:N-1);                     % Approximation weights for the dynamics
    R = diag([1 R]);                        % Approximation weights for the dynamics

    % Chebyshev polynomials difference
    S = zeros(N,N);
    S(1,1) = 1; 
    S(1,2) = -1/2;
    S(1,N) = (-1)^(order+1)/(order-1);
    for i = 2:order-1
        S(1,i+1) = (-1)^(i+1)*(1/(i-1)-1/(i+1));
        S(i,i-1) = 1; 
        S(i,i+1) = -1;
    end
    S(order,end-2) = 1;
    S(order,end) = -1;
    S(N,end-1) = 1;

    chi = zeros(N,size(x0,2));              % Boundary conditions term

    Ca = R*S*T.'*V;                         % Dynamics approximation matrix
    Cx = T*W;                               % State approximation matrix

    % Main loop
    while (GoOn && iter < maxIter)
        % Evaluate the dynamics 
        g = dynamics(tau.',x0); 

        % Compute the beta coefficients
        chi(1,:) = 2*x0(1,:);
        beta = Ca*g+chi; 

        % Compute the new state trajectory
        x = Cx*beta;

        % Convergence check 
        dX = x-x0; 
        dX = sqrt(dot(dX,dX,2));
        if (error < tol && norm(dX) < tol)
            GoOn = false;
        else
            iter = iter+1;
            x0 = x; 
            error = norm(dX);
        end
    end

    % Final output 
    state.State = ~GoOn; 
    state.Iterations = iter; 
    state.Error = error; 
end
