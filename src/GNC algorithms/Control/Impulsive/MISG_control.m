%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 07/10/22
% File: MISG_control.m 
% Issue: 0 
% Validated: 07/10/22

%% Multi impulsive, staging control %%
% This script contains the function to compute the control law by means of an MISG controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar Ln, the index of the libration point of interest
%         - scalar TOF, the time of flight for the rendezvous condition
%         - vector s0, initial conditions of both the target and the
%           relative particle
%         - string method, the solver to be used
%         - string integrator, the integrator to be used for the
%           variational equations
%         - scalar N, the number of impulses in the sequence
%         - scalar tol, the differential corrector scheme absolute
%           tolerance

% Output: - vector tspan, the sequence times
%         - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [tspan, Sc, dV, state] = MISG_control(mu, Ln, TOF, s0, method, integrator, N, tol)
    % Constants 
    m = 6;                      % Phase space dimension
    B = [zeros(3); eye(3)];     % Control input matrix

    stm_computation = integrator;

    switch (integrator)
        case 'RLLM'
            L = libration_points(mu);                       % System libration points
            gamma = L(end,Ln);                              % Characteristic distance of the libration point
            cn = legendre_coefficients(mu, Ln, gamma, 2);   % Legendre coefficient c_2 (equivalent to mu)
            c2 = cn(2);                                     % Legendre coefficient c_2 (equivalent to mu)
            Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];        % Potential linear term
            Omega = [0 2 0;-2 0 0; 0 0 0];                  % Coriolis force
            A = [zeros(m/2) eye(m/2); Sigma Omega];         % Jacobian matrix of the dynamics

        otherwise
            A = zeros(m);
    end
    
    % Sanity check on initial conditions dimension 
    if (size(s0,1) ~= 1)
        s0 = s0.';
    end
        
    % Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 

    tspan = linspace(0,TOF,N+1);                       % Integration time span
    dt = tspan(2)-tspan(1);                            % Time step
            
    % Preallocation of the impulses     
    dV = zeros(3*length(tspan),1);                     % Sequence of impulses
    sol = dV; 

    % Initial guess (null control)
    s0 = [s0 reshape(eye(m), [1 m^2])];
    [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);

    switch (method)
        case 'Numerical'
            % Differential corrector set up
            maxIter = 100;                                    % Maximum number of iterations
            GoOn(1) = true;                                   % Convergence boolean 
            iter(1) = 1;                                      % Initial iteration 

            % Outer loop
            while (GoOn(1) && iter(1) < maxIter) 
                % Inner loop constants 
                M = reshape(S(end,2*m+1:end), [m m]);         % Final STM
                epsilon = -S(end,m+1:2*m).';                  % Rendezvous error
                        
                % Inner loop (numerical solver) 
                sol = fsolve(@(x)optimal_sequence(S,M,B,epsilon,0.1,x), sol, optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt'));
                
                % Re-integrate the trajectory
                for i = 1:size(S,1)               
                    % New trajectory
                    U = [zeros(1,m+m/2) sol(1+m/2*(i-1):m/2*i).' zeros(1,m^2)];
                    [~, s] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), [0 dt], S(i,:)+U, options); 
        
                    % Update the state trajectory
                    S(i+1,:) = s(end,:);
                end
                S = S(1:end-1,:);
        
                % Convergence analysis 
                if (norm(dV-sol(1:m/2*size(S,1))) < tol)
                    GoOn(1) = false;                          % Stop the method
                else
                    dV = sol(1:m/2*size(S,1));                % Converged impulse sequence
                    iter(1) = iter(1)+1;                      % Update the iterations
                end
            end

            % Output
            dV = reshape(sol, [3 length(sol)/3]);             % Converged final impulses
            Sc = S;                                           % Control trajectory 
            state.State = ~GoOn(1);                           % Convergence boolean
            state.Iterations = iter;                          % Number of required iterations
            state.Error = norm(epsilon);                      % Final error

        case 'Sequential'
            % Differential corrector set up
            maxIter = 100;                                    % Maximum number of iterations
            GoOn(1) = true;                                   % Convergence boolean 
            iter(1) = 1;                                      % Initial iteration 

            % Outer loop
            while (GoOn(1) && iter(1) < maxIter) 
                % Inner loop constants 
                M = reshape(S(end,2*m+1:end), [m m]);         % Final STM
                epsilon = -S(end,m+1:2*m).';                  % Rendezvous error
                        
                % Inner loop (numerical solver) 
                if (iter(1) == 1)
                    sol = fsolve(@(x)optimal_sequence(S,M,B,epsilon,0.1,x), sol, optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt'));
                end

                GoOn(2) = true; 
                iter(2) = 1; 

                while (GoOn(2) && iter(2) < maxIter)
                    % Sequential update 
                    seq_eq = @(A,V0,x) (norm(V0)*x-A*norm(x)*V0);
                    for i = 1:length(tspan)-1
                        Phi0 = (M*reshape(S(i,2*m+1:end), [m m])^(-1)*B).';
                        Phi1 = (M*reshape(S(i+1,2*m+1:end), [m m])^(-1)*B).'; 
                        A = Phi1*pinv(Phi0);
                        sol(1+3*i:3*(i+1)) = fsolve(@(x)seq_eq(A, sol(1+3*(i-1):3*i), x), sol(1+3*i:3*(i+1)), optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt'));
                    end

                    % Initial impulse update
                    STM = zeros(m,m/2);             % Compute the complete STM
                    for i = 1:length(tspan)
                        % Initial STM
                        Phi0 = (M*reshape(S(1,2*m+1:end), [m m])^(-1)*B).';
                        Phi1 = (M*reshape(S(2,2*m+1:end), [m m])^(-1)*B).'; 
                        Prod = Phi1*pinv(Phi0);

                        % Cumulative STM
                        for j = 2:i-1
                            Phi0 = (M*reshape(S(j,2*m+1:end), [m m])^(-1)*B).';
                            Phi1 = (M*reshape(S(j+1,2*m+1:end), [m m])^(-1)*B).'; 
                            Prod = Phi1*pinv(Phi0)*Prod;
                        end 
                        Phi = M*reshape(S(i,2*m+1:end), [m m])^(-1);
                        stm = norm(sol(1+3*(i-1):3*i))*Phi*B*Prod;
                        STM = STM+stm;
                    end
                    e = STM*sol(1:3)/norm(sol(1:3))-epsilon;
                    G = STM*(eye(3)-sol(1:3)*sol(1:3).')/norm(sol(1:3))^3;

                    % Newton update 
                    ds = -pinv(G)*e;
                    sol(1:3) = sol(1:3)+ds;

                    % Convergence analysis 
                    if (norm(ds) < tol(2))
                        GoOn(2) = false;
                    else
                        iter(2) = iter(2)+1;
                    end
                end
                
                % Re-integrate the trajectory
                for i = 1:size(S,1)               
                    % New trajectory
                    U = [zeros(1,m+m/2) sol(1+m/2*(i-1):m/2*i).' zeros(1,m^2)];
                    [~, s] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), [0 dt], S(i,:)+U, options); 
        
                    % Update the state trajectory
                    S(i+1,:) = s(end,:);
                end
                S = S(1:end-1,:);
        
                % Convergence analysis 
                if (norm(epsilon) < tol(1))
                    GoOn(1) = false;                          % Stop the method
                else
                    iter(1) = iter(1)+1;                      % Update the iterations
                end
            end

            sol(end-2:end) = -s(end,10:12).';                 % Final impulse 
            S(end,10:12) = zeros(1,3);

            % Output
            dV = reshape(sol, [3 length(sol)/3]);             % Converged final impulses
            Sc = S;                                           % Control trajectory 
            state.State = ~GoOn(1);                           % Convergence boolean
            state.Iterations = iter;                          % Number of required iterations
            state.Error = norm(S(end,7:9));                   % Final error

        case 'Primal'
            s = S;                                              % Initial guess
            for i = 1:length(tspan)-1
                % Constants of the loop 
                M = reshape(s(end,2*m+1:end), [m m]);           % Final STM
                epsilon = -s(end,m+1:2*m).';                    % Rendezvous error
        
                % Solve the optimal problem
                solve_options = optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');
                aux = fsolve(@(x)optimal_sequence(s(i:end,2*m+1:end), M, B, epsilon, 0.1, x), sol(1+3*(i-1):end), solve_options);

                % Pruning
                dv = reshape(aux, 3, []);
                [dv, ~] = ISP_control(s(i:end,2*m+1:end), B, dv);
                sol(1+3*(i-1):3*i) = dv(1:3);
        
                % Update the propagation 
                switch (stm_computation)
                    case 'Numerical'
                        U = [zeros(1,m+m/2) sol(1+3*(i-1):3*i).' zeros(1,m^2)];
                        [~, aux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan(i:end), s(i,:)+U, options); 

                    case 'RLLM'
                        % Trajectory propagation
                        U = [zeros(1,m+m/2) sol(1+3*(i-1):3*i).'];
                        [~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan(i:end), s(i,1:2*m)+U, options); 

                        % STM propagation 
                        aux = zeros(length(tspan(i:end)), 2*m+m^2);
                        for j = 1:length(tspan(i:end))
                            STM = expm(A*tspan(i+(j-1)));
                            aux(j,:) = [S(j,:) reshape(STM, 1, m^2)];
                        end
                end

                if (length(tspan(i:end)) == 2)
                    s(i:end,:) = aux([1 end], :);
                else
                    s(i:end,:) = aux;
                end
            end  

            sol(end-2:end) = -s(end,10:12).';                 % Final impulse 
            s(end,10:12) = zeros(1,3);

            % Output
            dV = reshape(sol, [3 length(sol)/3]);             % Converged final impulses
            Sc = s;                                           % Control trajectory 
            state.State = true;                               % Convergence boolean
            state.Iterations = length(tspan);                 % Number of required iterations
            state.Error = norm(s(end,7:12));                  % Final error

        case 'Dual'
            s = S;                                              % Initial guess
            for i = 1:length(tspan)-1
                % Constants of the loop 
                M = reshape(s(end,2*m+1:end), [m m]);           % Final STM
                epsilon = -s(end,m+1:2*m).';                    % Rendezvous error
        
                % Solve the optimal problem
                [dv, ~, ~] = ADMM_sequence(s(i:end,2*m+1:end), M, B, epsilon, 'L1', 1e2, 1e3, 1e-4);

                % Pruning
                dv = reshape(dv, 3, []);
                [dv, ~] = ISP_control(s(i:end,2*m+1:end), B, dv);
                sol(1+3*(i-1):3*i) = dv(1:3);
        
                % Update the propagation 
                switch (stm_computation)
                    case 'Numerical'
                        U = [zeros(1,m+m/2) sol(1+3*(i-1):3*i).' zeros(1,m^2)];
                        [~, aux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan(i:end), s(i,:)+U, options); 

                    case 'RLLM'
                        % Trajectory propagation
                        U = [zeros(1,m+m/2) sol(1+3*(i-1):3*i).'];
                        [~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan(i:end), s(i,1:2*m)+U, options); 

                        % STM propagation 
                        aux = zeros(length(tspan(i:end)), 2*m+m^2);
                        for j = 1:length(tspan(i:end))
                            STM = expm(A*tspan(i+(j-1)));
                            aux(j,:) = [S(j,:) reshape(STM, 1, m^2)];
                        end
                end

                if (length(tspan(i:end)) == 2)
                    s(i:end,:) = aux([1 end], :);
                else
                    s(i:end,:) = aux;
                end
            end  

            sol(end-2:end) = -s(end,10:12).';                 % Final impulse 
            s(end,10:12) = zeros(1,3);

            % Output
            dV = reshape(sol, [3 length(sol)/3]);             % Converged final impulses
            Sc = s;                                           % Control trajectory 
            state.State = true;                               % Convergence boolean
            state.Iterations = length(tspan);                 % Number of required iterations
            state.Error = norm(s(end,7:12));                  % Final error

        otherwise
            error('No valid method was chosen');
    end
end

%% Auxiliary functions 
% Recursive optimal sequence
function [dV] = optimal_sequence(S, M, B, epsilon, dVmax, x)
    % State variables
    V = reshape(x(1:3*size(S,1)), [3 size(S,1)]);           % Optimal sequence
    dV = zeros(3*(size(S,1)-1),1);                          % Thrusting direction
    % t = x(3*size(S,1)+1:end);                             % Slack variable

    % Constants         
    m = size(B,1);          % Phase space dimension

    % Optimal thrusting direction recursion
    for i = 1:size(S,1)-1
        Phi0 = M*reshape(S(i,:), [m m])^(-1)*B;
        Phi1 = M*reshape(S(i+1,:), [m m])^(-1)*B;
        dV(1+3*(i-1):3*i,1) = norm(V(:,i))*V(:,i+1)-norm(V(:,i+1))*Phi1.'*pinv(Phi0.')*V(:,i);
    end

    % Rendezvous constraint
    Aeq = zeros(6,3*size(S,1));
    for i = 1:size(S,1)
        Aeq(1:m,1+m/2*(i-1):m/2*i) = M*reshape(S(i,:), [m m])^(-1)*B;
    end
    dV = [dV; Aeq*x(1:3*size(S,1))-epsilon];

    % Maximum impulse inequality
    % dV = [dV; dVmax^2-dot(V,V,1).'-t.^2];
end

% ADMM impulse sequence algorithm 
function [dV, lambda, cost] = ADMM_sequence(S, M, B, epsilon, cost_norm, maxIter, rho, tol)
    % Constants
    m = size(B,1);                            % Phase space dimension 

    % Preallocation and start up
    dV = zeros(3*(size(S,1)-1),maxIter);      % True impulse sequence 
    y = zeros(3*(size(S,1)-1),maxIter);       % Virtual L22 impulse sequence 
    mu = zeros(3*(size(S,1)-1),maxIter);      % Lagrange multiplier of the ADMM problem 
    GoOn = true;                              % Convergence boolean
    iter = 1;                                 % Initial iterations 
    J = zeros(1,maxIter);                     % Cost fuction for each iteration

    % ADMM optimization
    while (GoOn && iter < maxIter)
        % Quadratic convex optimization 
        epsilon_m = epsilon;
        for i = 1:(size(S,1)-1)
            Phi = M*reshape(S(i,:), [m m])^(-1)*B;
            epsilon_m = epsilon_m-Phi*rho/(2+rho)*(mu(1+3*(i-1):3*i,iter)+dV(1+3*(i-1):3*i,iter));
        end
        [y(:,iter+1), lambda, ~] = lagrange_sequence(S, M, B, epsilon_m);

        % Feedforward term
        y(:,iter+1) = y(:,iter+1)+rho/(2+rho)*(mu(:,iter)+dV(:,iter));
        dVy = reshape(y, 3, []); 
        cost = sum(dot(dVy,dVy,1));

        % Proximal operator 
        V = zeros(1,(size(S,1)-1));
        for i = 1:(size(S,1)-1)
            dV(1+3*(i-1):3*i,iter+1) = proximal_operator(y(1+3*(i-1):3*i,iter+1)-mu(1+3*(i-1):3*i,iter), 1/rho, cost_norm);
            V(i) = norm(dV(1+3*(i-1):3*i,iter+1));
        end

        % Lagrange multiplier update
        mu(:,iter+1) = mu(:,iter) + dV(:,iter)-y(:,iter);

        % Total cost 
        res = reshape(dV(:,iter+1)-y(:,iter+1), [3 (size(S,1)-1)]);
        J(iter+1) = sum(V)+cost+rho/2*sum(sqrt(dot(res,res,1)))+dot(mu(:,iter+1),dV(:,iter+1)-y(:,iter+1));
        cost = sum(V);

        % Convergence analysis 
        if (abs(J(iter+1)-J(iter)) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end

    dV = dV(:,iter);
end

% Lagrange optimal equation 
function [dV, lambda, cost] = lagrange_sequence(S, M, B, epsilon)
    % Rendezvous constraint
    m = 6;
    Aeq = zeros(size(M));
    for i = 1:size(S,1)
        C = M*reshape(S(i,:), [m m])^(-1)*B;
        Aeq = Aeq+C*C.';
    end

    % Preallocation
    dV = zeros(3*(size(S,1)-1),1);      % True impulse sequence

    % Solve for the Lagrange multiplier
    lambda = Aeq\epsilon;

    % Impulse sequence 
    for i = 1:size(S,1)-1
        Phi = M*reshape(S(i,:), [m m])^(-1)*B;
        dV(1+3*(i-1):3*i,1) = Phi.'*lambda;
    end

    % Final cost 
    cost = 0;
    for i = 1:size(S,1)-1
        cost = cost + norm(dV(1+3*(i-1):3*i,1));
    end
end