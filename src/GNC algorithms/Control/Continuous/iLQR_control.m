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

function [tspan, Sc, u, state] = iLQR_control(mu, tf, s0, GNC)
    % Constants of the model 
    n = 6;                        % Dimension of the state vector
    tol = 1e-10;                  % Convergence tolerance

    b = [zeros(3); eye(3)];       % Control input matrix
    A = [zeros(n,n)];             % PID-structured Jacobian matrix

    % Controller definition
    Q = GNC.Control.LQR.Q;        % Penalty matrix on the state error
    R = GNC.Control.LQR.M;        % Penalty matrix on the control effort

    % Generate the initial LQR/SDRE guess 
    s0 = [s0 reshape(eye(n), [1 n^2])];                                                                   % Variational initial conditions
    dt = 1e-3;                                                                                            % Nondimensional time step
    tspan = 0:dt:tf;                                                                                      % Integration time span
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);                                                % Integration tolerances
    [~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);          % Initial guess for the state trajectory
    
    Q = Q(1:n,1:n);
    u = zeros(3,size(Sc,1));

    V = 0; 
    for i = 1:length(tspan)
        V = V+0.5*Sc(i,n+1:2*n)*Q*Sc(i,n+1:2*n).'+0.5*u(:,i).'*R*u(:,i);
    end

    % Preallocation 
    du = zeros(size(b,2),size(Sc,1));                   % Update to the control law 
    ds = zeros(n,size(Sc,1));                           % Update to the rendezvous trajectory
    v = zeros(n,size(Sc,1));                            % Feedforward term to the rendezvous trajectory

    K = zeros(size(u,1),size(ds,1)*size(Sc,1));         % LQR gain
    Ku = zeros(size(u,1),size(u,1)*size(Sc,1));         % Feedback gain
    Kv = zeros(size(u,1),size(v,1)*size(Sc,1));         % Feeforward gain
    
    % iLQR iterative loop 
    GoOn = true;                            % Convergence boolean 
    maxIter = 100;                          % Maximum number of iterations
    iter = 1;                               % Initial iteration 

    % Boundary conditions 
    ds(:,1) = zeros(n,1);                   % Initial update

    while (GoOn && iter < maxIter)
        % Boundary conditions
        S = Q;                                  % Penalty weight matrix final conditions
        v(:,end) = Q*Sc(end,n+1:2*n).';         % Final feedforward term

        for i = size(Sc,1)-1:-1:1
            % Linearization 
            if (i > 1)
                Phi0 = reshape(Sc(i-1,2*n+1:end), [n n]);   % STM from 0 to k-1
            else
                Phi0 = eye(n);                              % STM from 0 to k-1
            end
            Phi = reshape(Sc(i,2*n+1:end), [n n]);          % STM from 0 to k
            A(1:n,1:n) = Phi*Phi0^(-1);                     % Relative dynamics linear model
            B = A*b;                                        % Control input matrix

            % Controller gains
            K(:,1+size(ds,1)*(i-1):size(ds,1)*i) = (B.'*S*B+R)^(-1)*B.'*S*A;
            Kv(:,1+size(v,1)*(i-1):size(v,1)*i) = (B.'*S*B+R)^(-1)*B.';
            Ku(:,1+size(u,1)*(i-1):size(u,1)*i) = (B.'*S*B+R)^(-1)*R;

            % Backward pass recursion
            S = A.'*S*(A-B*K(:,1+size(ds,1)*(i-1):size(ds,1)*i))+Q; 
            v(:,i) = (A-B*K(:,1+size(ds,1)*(i-1):size(ds,1)*i)).'*v(:,i+1)-K(:,1+size(ds,1)*(i-1):size(ds,1)*i).'*R*u(:,i)+Q*Sc(i,n+1:2*n).';
        end

        for i = 1:size(Sc,1)-1
            % Linearization 
            if (i > 1)
                Phi0 = reshape(Sc(i-1,2*n+1:end), [n n]);   % STM from 0 to k-1
            else
                Phi0 = eye(n);                              % STM from 0 to k-1
            end
            Phi = reshape(Sc(i,2*n+1:end), [n n]);          % STM from 0 to k
            A(1:n,1:n) = Phi*Phi0^(-1);                     % Relative dynamics linear model
            B = A*b;                                        % Control input matrix

            % Backward pass
            du(:,i) = -(K(:,1+size(ds,1)*(i-1):size(ds,1)*i)*ds(:,i)+Ku(:,1+size(u,1)*(i-1):size(u,1)*i)*u(:,i)+Kv(:,1+size(v,1)*(i-1):size(v,1)*i)*v(:,i+1));
            ds(:,i+1) = A*ds(:,i)+B*du(:,i);
        end

        % Forward pass        
        u = u+du;
%         dynamics = @(t,s)(vnlr_model(mu, u, t.', s.').');
%         [tau, Sc, state] = MCPI([tspan(1) tspan(end)], tau, Sc, dynamics, N, 1e-12);
        for i = 1:size(Sc,1)-1
            [~,s] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), [0 dt], Sc(i,:)+[zeros(1,9) u(:,i).' zeros(1, n^2)], options);
            Sc(i+1,:) = s(end,:);
        end

        Vp = 0; 
        for i = 1:length(tspan)
            Vp = Vp+0.5*Sc(i,n+1:2*n)*Q*Sc(i,n+1:2*n).'+0.5*u(:,i).'*R*u(:,i);
        end

        % Convergence check 
        dV = abs(Vp-V);
        if (dV < tol) 
            GoOn = false;
        else
            iter = iter+1; 
            V = Vp;
        end
    end

    % Final output
    state.State = ~GoOn;
    state.Iterations = iter; 
    state.Error = dV;
end

function [Sc, u, state] = iLQR_LT(mu, tf, s0, GNC)
    % Constants of the model 
    n = 6;                        % Dimension of the state vector
    tol = 1e-5;                   % Convergence tolerance

    B = [zeros(3); eye(3)];               % Control input matrix
    A = [zeros(n,n)];          % PID-structured Jacobian matrix

    % Controller definition
    Q = GNC.Control.LQR.Q;                          % Penalty matrix on the state error
    R = GNC.Control.LQR.M;                          % Penalty matrix on the control effort

    % Generate the initial LQR/SDRE guess 
%     dt = 1e-3;                                                                                            % Nondimensional time step
%     tspan = 0:dt:tf;                                                                                      % Integration time span
%     options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);                                                % Integration tolerances
%     [t, Sc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);    % Initial guess for the state trajectory
%     [~, ~, u] = GNC_handler(GNC, Sc(:,1:n), Sc(:,n+1:end), t);                                            % Initial guess for the control vector

    N = 1e3;
    tau = flip(cos((0:N)*pi/N));
    Sc = zeros(N+1,2*n);
    Sc(1,:) = s0(1:2*n); 
    u = zeros(3,N+1);

    dynamics = @(t,s)(vnlr_model(mu, u, t.', s.').');
    [tspan, Sc, state] = MCPI([0 tf], tau, Sc, dynamics, N, 1e-12);

    Sc = Sc(:,1:2*n);
    Q = Q(1:n,1:n);

    V = 0; 
    for i = 1:length(tspan)
        V = V+0.5*Sc(i,n+1:end)*Q*Sc(i,n+1:end).'+0.5*u(:,i).'*R*u(:,i);
    end

    % Preallocation 
    du = zeros(size(B,2),size(Sc,1));                   % Update to the control law 
    ds = zeros(n,size(Sc,1));                         % Update to the rendezvous trajectory
    v = zeros(n,size(Sc,1));                          % Feedforward term to the rendezvous trajectory

    K = zeros(size(u,1),size(ds,1)*size(Sc,1));         % LQR gain
    Ku = zeros(size(u,1),size(u,1)*size(Sc,1));         % Feedback gain
    Kv = zeros(size(u,1),size(v,1)*size(Sc,1));         % Feeforward gain
    
    % iLQR iterative loop 
    GoOn = true;                            % Convergence boolean 
    maxIter = 100;                          % Maximum number of iterations
    iter = 1;                               % Initial iteration 

    % Boundary conditions 
    ds(:,1) = zeros(n,1);                   % Initial update

    while (GoOn && iter < maxIter)
        % Boundary conditions
        S = Q;                                  % Penalty weight matrix final conditions
        v(:,end) = Q*Sc(end,n+1:end).';         % Final feedforward term

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
            du(:,i) = -(K(:,1+size(ds,1)*(i-1):size(ds,1)*i)*ds(:,i)+Ku(:,1+size(u,1)*(i-1):size(u,1)*i)*u(:,i)+Kv(:,1+size(v,1)*(i-1):size(v,1)*i)*v(:,i+1));
            ds(:,i+1) = A*ds(:,i)+B*du(:,i);
        end

        % Forward pass        
        u = u+du;
        dynamics = @(t,s)(vnlr_model(mu, u, t.', s.').');
        [tau, Sc, state] = MCPI([tspan(1) tspan(end)], tau, Sc, dynamics, N, 1e-12);

        Vp = 0; 
        for i = 1:length(tspan)
            Vp = Vp+0.5*Sc(i,n+1:end)*Q*Sc(i,n+1:end).'+0.5*u(:,i).'*R*u(:,i);
        end

        % Convergence check 
        dV = abs(Vp-V);
        if (dV < tol) 
            GoOn = false;
        else
            iter = iter+1; 
            V = Vp;
        end
    end

    % Final output 
    state.State = ~GoOn; 
    state.Iterations = iter; 
    state.Error = ds; 
end