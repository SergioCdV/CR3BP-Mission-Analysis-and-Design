%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/09/22
% File: iLQR_control.m 
% Issue: 0 
% Validated: 24/09/22

%% Iterative Linear Quadratic Regulator Control %%
% This script contains the function to compute the control law by means of an iLQR controller.

% Inputs: - string model, selecting the linear model to compute the linear state
%           transition matrix of the system
%         - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - array St, the target orbit reference state 
%         - array Sg, the guidance law to follow
%         - array Sn, the system state
%         - scalar Ln, the libration point number. It may be left 0 if the
%           target orbit is not librating around any Lagrange point
%         - scalar gamma, the relative distance of the Ln point to the
%           nearest primary. Again, it may be left as 0 if needed
%         - matrices Q and M, penalizing on the state error and the control
%           effort

% Output: - array Sc, the rendezvous relative trajectory
%         - array u, containing  column-wise the required control law
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [Sc, u, state] = iLQR_control(mu, tf, s0, GNC)
    % Constants of the model 
    n = 6;                        % Dimension of the state vector
    tol = 1e-5;                   % Convergence tolerance

    B = [zeros(3); eye(3); zeros(3)];               % Control input matrix
    A = [zeros(n,n+3); eye(3) zeros(3,n)];          % PID-structured Jacobian matrix

    % Controller definition
    Q = GNC.Control.LQR.Q;                          % Penalty matrix on the state error
    R = GNC.Control.LQR.M;                          % Penalty matrix on the control effort

    % Generate the initial LQR/SDRE guess 
    dt = 1e-3;                                                                                            % Nondimensional time step
    tspan = 0:dt:tf;                                                                                      % Integration time span
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);                                                % Integration tolerances
    [t, Sc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);    % Initial guess for the state trajectory
    [~, ~, u] = GNC_handler(GNC, Sc(:,1:n), Sc(:,n+1:end), t);                                            % Initial guess for the control vector

    % Preallocation 
    du = zeros(3,size(Sc,1));                           % Update to the control law 
    ds = zeros(n+3,size(Sc,1));                         % Update to the rendezvous trajectory
    v = zeros(n+3,size(Sc,1));                          % Feedforward term to the rendezvous trajectory

    K = zeros(size(u,1),size(ds,1)*size(Sc,1));         % LQR gain
    Ku = zeros(size(u,1),size(u,1)*size(Sc,1));         % Feedback gain
    Kv = zeros(size(u,1),size(v,1)*size(Sc,1));         % Feeforward gain
    
    % iLQR iterative loop 
    GoOn = true;                            % Convergence boolean 
    maxIter = 100;                          % Maximum number of iterations
    iter = 1;                               % Initial iteration 

    % Boundary conditions 
    ds(:,1) = zeros(n+3,1);                 % Initial update
    S = Q;                                  % Penalty weight matrix final conditions
    v(:,end) = Q*Sc(end,n+1:end).';         % Final feedforward term

    while (GoOn && iter < maxIter)
        for i = size(Sc,1)-1:-1:1
            % Linearization 
            A(1:n,1:n) = rel_jacobian(mu, Sc(i,:).');         % Jacobian of the relative dynamics

            % Controller gains
            K(:,1+size(ds,1)*(i-1):size(ds,1)*i) = (B.'*S*B+R)^(-1)*B.'*S*A;
            Kv(:,1+size(v,1)*(i-1):size(v,1)*i) = (B.'*S*B+R)^(-1)*B.';
            Ku(:,1+size(u,1)*(i-1):size(u,1)*i) = (B.'*S*B+R)^(-1)*R;

            % Backward pass recursion
            S = A.'*S*(A-B*K(:,1+size(ds,1)*(i-1):size(ds,1)*i))+Q; 
            v(:,i) = (A-B*K(:,1+size(ds,1)*(i-1):size(ds,1)*i)).'*v(:,i+1)-K(:,1+size(ds,1)*(i-1):size(ds,1)*i).'*R*u(:,i)+Q*Sc(i,n+1:end).';
        end

        for i = 1:size(Sc,1)-1
            % Linearization 
            A(1:n,1:n) = rel_jacobian(mu, Sc(i,:).');         % Jacobian of the relative dynamics 

            % Backward pass
            du(:,i) = -(K(:,1+size(ds,1)*(i-1):size(ds,1)*i)*ds(:,i)+Ku(:,1+size(u,1)*(i-1):size(u,1)*i)*u(:,i)+Kv(:,1+size(v,1)*(i-1):size(v,1)*i)*v(:,i));
            ds(:,i+1) = A*ds(:,i)+B*du(:,i);
        end

        % Forward pass        
        u = u+du;
        Sc(:,n+1:end) = Sc(:,n+1:end)+ds.';
%         for i = 1:size(Sc,1)-1
%             F(:,1) = [nlr_model(mu, true, false, false, 'Encke', tspan(i), Sc(i,1:2*n).'); Sc(i,n+1:n+3).']+[zeros(9,1); u(:,i); zeros(3,1)];
%             aux = [Sc(i+1,1:n).'; Sc(i+1,n+1:end).'+F(n+1:end,1)*dt];
%             F(:,2) = [nlr_model(mu, true, false, false, 'Encke', tspan(i+1), aux(1:2*n)); aux(n+1:n+3)]+[zeros(9,1); u(:,i+1); zeros(3,1)];
%             Sc(i+1,n+1:end) = Sc(i,n+1:end) + dt*sum(F(n+1:end,:),2).'/2;
%         end

        % Convergence check 
        dS = sqrt(dot(du,du,1)); 
        max(dS)
        if (max(dS) < tol) 
            GoOn = false;
        else
            iter = iter+1; 
        end
    end

    % Final output 
    state.State = ~GoOn; 
    state.Iterations = iter; 
    state.Error = ds; 
end