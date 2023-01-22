

function [tspan, St, u] = PFKS_def(mu, TOF, T, initial, alpha_min, GNC)
    % Constants 
    n = 6; 
    Tmax = GNC.Tmax; 

    % Guidance 
    Jref = jacobi_constant(mu, initial(1,1:n).');            % Guidance law for the Jacobi constant

    tspan = 0:1e-3:T;                                        % Maneuver time span
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    [~, St] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, [initial(1:n) reshape(eye(n), 1, [])], options);

    [F0, lambda] = eig(reshape(St(end,n+1:end), [n n]));
    J = diag(log(diag(lambda))/T);
    for i = 1:size(F0,2)
        F0(:,i) = F0(:,i)/lambda(i,i);
    end
  
    tspan = linspace(0, TOF, 60);
    F = zeros(6,6*length(tspan)); 
    for i = 1:length(tspan)
        F(:,1+n*(i-1):n*i) = reshape(St(i,n+1:end), [n n])*F0*expm(-J*mod(tspan(i),T));
    end

    [~, St] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, [initial(1:n) reshape(eye(n), 1, [])], options);

    % New initial conditions 
    initial = [F(:,1:6)^(-1)*initial(n+1:end).'; jacobi_constant(mu, St(1,1:n).'+initial(n+1:end).')];
  
    % Domain definition
    m = length(tspan)-1;
    [tau, ~, ~, D] = quadrature(3, m, 'Legendre');

    % Initial guess
    initial = real(initial([1 2 end]));
    P0 = repmat(initial, 1, length(tau));
    x0 = reshape(P0, [], 1);
    x0 = [x0; zeros(length(tau),1); 2; 2];
    
    % Upper and lower bounds 
    P_lb = [-Inf*ones(size(reshape(P0,[],1))); -Inf*ones(length(tau),1); -Inf; -Inf];
    P_ub = [Inf*ones(size(reshape(P0,[],1))); Inf*ones(length(tau),1); Inf; Inf];
    
    % Objective function
    objective = @(x)cost_function(tau, x);

    % Non-linear constraints
    cons = @(x)nonlcon(mu, initial, Jref, alpha_min, St, F, J, D, tspan, tau, x);
    
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, alpha_2, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, cons, options);
    
    % Solution 
    C = reshape(sol(1:3*(m+1)), [size(P0,1) size(P0,2)]);             % Optimal control points
    lambda = reshape(sol(3*length(tau)+1:end-2), [1 size(P0,2)]);   % Optimal control points
    lambda0 = sol(end-1:end);

    Lambda = zeros(3,length(tau));
    u = zeros(3,length(tau)); 

    for i = 1:length(tau)
        dJ = jacobi_gradient(mu, St(i,1:6).');
        Lambda(:,i) = [expm(-J(1:2,1:2).'*tspan(i))*lambda0; lambda(i)];

        E = F(:,1+6*(i-1):6*i)^(-1);
        G = real([E(1:2,4:6); dJ(4:6).']);  % Control operator

        u(:,i) = -Tmax*G.'*Lambda(:,i)/norm(G.'*Lambda(:,i));
    end

    fprintf('Exit flag: %i\n', exitflag)
    if (exitflag ~= 1)
        fprintf("Exit messsage: %s", output.message);
    end

    fprintf("Number of iterations: %i\n", output.iterations);
    fprintf("Number of function evaluations: %i\n", output.funcCount);
    fprintf("Constraint violation: %f \n", output.constrviolation);
    
    % Results 
    figure
    plot(tau, C); 
    legend('$\alpha_u$', '$\alpha_s$', '$J$')
    xlabel("$t$")
    ylabel("$\mathbf{s}$")
    grid on;

    figure
    hold on
    plot(tau, u)
    xlabel("$t$")
    ylabel("$\mathbf{u}$")
    grid on;
end

%% Auxiliary functions
% Control cost 
function [cost] = cost_function(tau, x)
    % Reshaping 
    C = reshape(x(1:3*length(tau)), 3, length(tau));

    % Mayer term 
    M = -C(2,end); 

    % Lagrange term 
    L = 0; 

    % Total cost 
    cost = M+L;
end

% Constraints
function [c, ceq] = nonlcon(mu, initial, Jref, alpha_min, St, F, J, D, t, tau, x)
    % Reshaping 
    C = reshape(x(1:3*length(tau)), 3, length(tau));                   % State variables
    Lambda_e = reshape(x(3*length(tau)+1:end-2), 1, length(tau));      % Lagrange multiplier of the energy equation
    Lambda_alpha = x(end-1:end);                                       % Initial conditions for the stationkeeping costates

    % Evolution of the costates 
    Lambda = zeros(3,length(t));
    for i = 1:length(t)
        Lambda(:,i) = [Lambda_alpha.*[exp(-J(1,1)*t(i)); exp(-J(2,2)*t(i))]; Lambda_e(i)];
    end

    % Boundary conditions constraints 
    ceq = [C(:,1)-initial; C(3,end)-Jref; Lambda(2,end)+C(2,end)];

    % Dynamic constraints 
    dyn(1:3,:) = C(1:3,:)*D.';
    dyn(4,:) = Lambda(3,:)*D.';
    u = zeros(3,length(tau));

    for i = 1:length(tau)
        % Matrices
        dJ = jacobi_gradient(mu, St(i,1:6).');

        K = dJ.'*F(:,1+6*(i-1):6*i)*J;               % Gradient of the Jacobi constant with respect to the relative state variables
        A = [J(1:2,1:2) zeros(2,1); K(1:2) 0];       % State dynamics

        E = F(:,1+6*(i-1):6*i)^(-1);
        G = real([E(1:2,4:6); dJ(4:6).']);  % Control operator

        u(:,i) = -G.'*Lambda(:,i);
        norm_u = norm(u(:,i)); 
        
        % Dynamics
        dyn(1:3,i) = dyn(1:3,i)-(t(end)/2)*(A*C(:,i)+G*u(:,i)/norm_u);
        aux = -A.'*Lambda(:,i);
        dyn(4,i) = dyn(4,i)-(t(end)/2)*aux(3);

        u(:,i) = u(:,i)/norm_u;
    end

    ceq = [ceq; max(abs(reshape(dyn, [], 1)))];
    
    % Inequality constraints 
    c = abs(C(1,end))-alpha_min;
end