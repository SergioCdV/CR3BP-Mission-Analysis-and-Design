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
%         - scalar TOF, the time of flight for the rendezvous condition
%         - vector s0, initial conditions of both the target and the
%           relative particle
%         - string method, the solver to be used
%         - scalar N, the number of impulses in the sequence
%         - scalar tol, the differential corrector scheme absolute
%           tolerance

% Output: - vector tspan, the sequence times
%         - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [tspan, Sc, dV, state] = MISG_control(mu, TOF, s0, method, N, tol)
    % Constants 
    m = 6;                      % Phase space dimension
    B = [zeros(3); eye(3)];     % Control input matrix
    
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
        case 'MPC'
            s = S;                                              % Initial guess
            for i = 1:length(tspan)-1
                % Constants of the loop 
                M = reshape(s(end,2*m+1:end), [m m]);           % Final STM
                epsilon = -s(end,m+1:2*m).';                    % Rendezvous error
        
                % Solve the optimal problem
                aux = fsolve(@(x)optimal_sequence(s(i:end,:),M,B,epsilon,0.1,x), sol(1+3*(i-1):end), optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt'));
                sol(1+3*(i-1):3*i) = aux(1:3);
        
                % Update the propagation 
                U = [zeros(1,m+m/2) sol(1+3*(i-1):3*i).' zeros(1,m^2)];
                [~, aux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan(i:end), s(i,:)+U, options); 
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
            state.Error = norm(epsilon);                      % Final error

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
                    if (norm(ds) < tol(1))
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
                if (norm(dV-sol(1:m/2*size(S,1))) < tol(2))
                    GoOn(1) = false;                          % Stop the method
                else
                    dV = sol(1:m/2*size(S,1));                % Converged impulse sequence
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
            state.Error = norm(epsilon);                      % Final error

        otherwise
            error('No valid method was chosen');
    end
end

%% Auxiliary functions 
function [dV] = optimal_sequence(S, M, B, epsilon, dVmax, x)
    % State variables
    V = reshape(x(1:3*size(S,1)), [3 size(S,1)]);           % Optimal sequence
    % t = x(3*size(S,1)+1:end);                             % Slack variable

    % Optimal thrusting direction
    dV = zeros(3*(size(S,1)-1),1);
    for i = 1:size(S,1)-1
        Phi0 = M*reshape(S(i,13:end), [6 6])^(-1)*B;
        Phi1 = M*reshape(S(i+1,13:end), [6 6])^(-1)*B;
        dV(1+3*(i-1):3*i,1) = norm(V(:,i))*V(:,i+1)-norm(V(:,i+1))*Phi1.'*pinv(Phi0.')*V(:,i);
    end

    % Rendezvous constraint
    m = 6;
    Aeq = zeros(6,3*size(S,1));
    M = reshape(S(end,2*m+1:end), [m m]);
    for i = 1:size(S,1)
        Aeq(1:m,1+m/2*(i-1):m/2*i) = M*reshape(S(i,2*m+1:end), [m m])^(-1)*B;
    end
    dV = [dV; Aeq*x(1:3*size(S,1))-epsilon];

    % Maximum impulse inequality
    % dV = [dV; dVmax^2-dot(V,V,1).'-t.^2];
end