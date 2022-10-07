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
%         - scalar tol, the differential corrector scheme absolute
%           tolerance

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [Sc, dV, state] = MISG_control(mu, TOF, s0, tol)
    % Constants 
    m = 6;                      % Phase space dimension
    B = [zeros(3); eye(3)];     % Control input matrix
    
    % Sanity check on initial conditions dimension 
    if (size(s0,1) ~= 1)
        s0 = s0.';
    end
        
    % Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    
    dt = 1e-2;                                         % Integration time step  
    tspan = 0:dt:TOF;                                  % Integration time span
    
    % Differential corrector set up
    maxIter = 100;                                     % Maximum number of iterations
    GoOn(1) = true;                                    % Convergence boolean 
    iter(1) = 1;                                       % Initial iteration 
        
    % Preallocation of the impulses 
    dV = zeros(length(tspan),1);                       % Vector of impulses
    u = zeros(3,length(tspan));                        % Thrusting direction
    error = zeros(length(tspan),1);                    % Constraint vector
    dVc = 10*ones(size(dV));                           % Converged impulse sequence

    % Initial guess based on the TISS solution
    [~, dv, ~] = TISS_control(mu, TOF, s0, tol, 'Position', true);
    dV = repmat(1e-1, length(tspan), 1);
    dV(end) = norm(dv(:,end));
              
    s0 = [s0 reshape(eye(m), [1 m^2])];
    [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0+[zeros(1,m+m/2) dv(:,1).' zeros(1,m^2)], options);
    S(end,10:12) = S(end,10:12)+dv(:,end).';
    
    % Outer loop
    while (GoOn(1) && iter(1) < maxIter) 
        % Inner loop constants 
        M = reshape(S(end,2*m+1:end), [m m]); 

        % Initial error 
        epsilon = -S(end,m+1:2*m).';

        % Re-initialization
        GoOn(2) = true; 
        iter(2) = 1; 

        % Inner loop (optimzing planning) 
        while (GoOn(2) && iter(2) < maxIter)
            % Compute the Lagrange multiplier 
            STM = zeros(m); 
            for i = 1:length(dV)
                Phi = M*reshape(S(i,2*m+1:end), [m m])^(-1);
                STM = STM + dV(i)*Phi*B*(B.'*Phi.');
            end
            lambda = exp(sum(dV))*STM\epsilon; 

            % Compute the thrusting directions 
            u(:,1) = B.'*M.'*lambda/exp(sum(dV));
            for i = 1:length(dV)-1
                Phi0 = (M*reshape(S(i,2*m+1:end), [m m])^(-1)*B).';
                Phi =  (M*reshape(S(i+1,2*m+1:end), [m m])^(-1)*B).';
                u(:,i+1) = dV(i+1)/dV(i)*Phi*pinv(Phi0)*u(:,i);
            end

            % Compute the error
            for i = 1:length(dV)
                Phi = M*reshape(S(i,2*m+1:end), [m m])^(-1)*B;
                error(i) = exp(sum(dV))-dot(u(:,i), Phi.'*lambda);
            end

            % Compute the Newton update
            ds = -error / exp(sum(dV)); 

            % Update the velocity norms 
            dV = dV + ds;

            % Convergence
            if (norm(ds) < tol)
                GoOn = false; 
            else
                iter(2) = iter(2) + 1; 
            end
        end
        
        % Re-integrate the trajectory
        for i = 1:size(dV,1)               
            % New trajectory
            U = [zeros(1,m+m/2) dV(i)*u(:,i).' zeros(1,m^2)];
            [~, s] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), [0 dt], S(i,:)+U, options); 

            % Update the state trajectory
            S(i+1,:) = s(end,:);
        end

        % Convergence analysis 
        if (norm(dVc-dV) < tol)
            GoOn(1) = false;                          % Stop the method
        else
            dVc = dV;                                 % Converged impulse sequence
            iter(1) = iter(1)+1;                      % Update the iterations
        end
    end
            
    % Output       
    dV = dV.*u;                                       % Converged final impulses
    Sc = S;                                           % Control trajectory 
    state.State = ~ (GoOn(1) && GoOn(2));             % Convergence boolean
    state.Iterations = max(iter);                     % Number of required iterations
    state.Error = norm(error);                        % Final error
end