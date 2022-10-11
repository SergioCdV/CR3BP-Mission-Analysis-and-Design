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

function [tspan, Sc, u, state] = iLQR_control(mu, tf, s0, GNC, tol)
    % Constants of the model 
    n = 6;                        % Dimension of the state vector

    % Controller definition
    mode = GNC.Control.iLQR.Mode; % Control options
    Q = GNC.Control.LQR.Q;        % Penalty matrix on the state error
    R = GNC.Control.LQR.M;        % Penalty matrix on the control effort

    m = size(Q,1);                % Dimension of the possibly augmented problem

    % Discrete dynamics definition 
    b = zeros(m,3);               % Control input matrix
    b(4:6,:) = eye(3);            % Control input matrix
    A = zeros(m);                 % PID-structured Jacobian matrix

    % Generate the initial LQR/SDRE guess 
    s0 = [s0 reshape(eye(m), [1 m^2])];                            % Variational initial conditions
    dt = 1e-2;                                                     % Nondimensional time step
    tspan = 0:dt:tf;                                               % Integration time span
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);         % Integration tolerances

    switch (mode)
        case 'Discrete'
            [~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);  % Initial guess for the state trajectory
            u = zeros(3,size(Sc,1));                                                                      % Initial guess for the control law

        case 'Continuous'
            N = length(tspan); 
            tau = flip(cos((0:N)*pi/N));
            tspan = (tspan(end)-tspan(1))/2*tau+(tspan(end)+tspan(1))/2; 
            [~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan, s0, options);  % Initial guess for the state trajectory
            [~, ~, u] = GNC_handler(GNC, Sc(:,1:n), Sc(:,n+1:n+m), tspan);                                     % Initial guess for the control law

        otherwise 
            error('No valid iLQR mode was selected');
    end

    % Initial cost function evaluation
    V = 0;
    switch (mode)
        case 'Discrete'
            for i = 1:length(tspan)
                V = V+0.5*Sc(i,n+1:n+m)*Q*Sc(i,n+1:n+m).'+0.5*u(:,i).'*R*u(:,i);
            end
        case 'Continuous'           
            for i = 1:length(tspan)
                V = V+0.5*Sc(i,n+1:n+m)*Q*Sc(i,n+1:n+m).'+0.5*dt*u(:,i).'*R*u(:,i);
            end
    end

    % Inequality penalty evaluations
    rho = 0; 
    lambda = zeros(1,size(u,2));
    C = sqrt(dot(u,u,1))-0.01;
    V = V+dot(lambda,C)+0.5*rho*dot(C,C);

    % Preallocation 
    du = zeros(size(b,2),size(Sc,1));                   % Update to the control law 
    ds = zeros(m,size(Sc,1));                           % Update to the rendezvous trajectory
    v = zeros(m,size(Sc,1));                            % Feedforward term to the rendezvous trajectory
    d = zeros(3,size(Sc,1));                            % Feedforward term to the rendezvous trajectory

    K = zeros(size(u,1),size(ds,1)*size(Sc,1));         % LQR gain
    Ku = zeros(size(u,1),size(u,1)*size(Sc,1));         % Feedback gain
    Kv = zeros(size(u,1),size(v,1)*size(Sc,1));         % Feeforward gain
    
    % iLQR iterative loop 
    GoOn(1) = true;                         % Outer loop convergence boolean 
    GoOn(2) = true;                         % Inner loop convergence boolean 
    maxIter = 100;                          % Maximum number of iterations
    iter(1) = 1;                            % Outer loop initial iteration 
    iter(2) = 1;                            % Inner loop initial iteration 

    % Boundary conditions 
    ds(:,1) = zeros(m,1);                   % Initial update

    while (GoOn(1) && iter(1) < maxIter)
        % Inner LQR loop
        GoOn(2) = true;                     % Convergence boolean reset
        iter(2) = 1;                        % Iterations reset

        while (GoOn(2) && iter(2) < maxIter)
            % Boundary conditions
            S = Q;                              % Penalty weight matrix final conditions
            v(:,end) = Q*Sc(end,n+1:n+m).';     % Final feedforward term
    
            for i = size(Sc,1)-1:-1:1
                % Linearization 
                if (i > 1)
                    Phi0 = reshape(Sc(i-1,n+m+1:end), [m m]);   % STM from 0 to k-1
                else
                    Phi0 = eye(m);                              % STM from 0 to k-1
                end
                Phi = reshape(Sc(i,n+m+1:end), [m m]);          % STM from 0 to k
                A(1:m,1:m) = Phi*Phi0^(-1);                     % Relative dynamics linear model
    
                switch (mode)
                    case 'Discrete'
                        B = A*b;                                                        % Control input matrix
                    case 'Continuous'
                        B = (dt/2)*(eye(m)+A)*b;                              % Control input matrix
                end
    
                % Controller gains
%                 K(:,1+size(ds,1)*(i-1):size(ds,1)*i) = (B.'*S*B+R)^(-1)*B.'*S*A;
%                 Kv(:,1+size(v,1)*(i-1):size(v,1)*i) = (B.'*S*B+R)^(-1)*B.';
%                 Ku(:,1+size(u,1)*(i-1):size(u,1)*i) = (B.'*S*B+R)^(-1)*R;

                lx = Q*Sc(i,n+1:n+m).';
                lu = R*u(:,i);
                lxx = Q;
                luu = R;
                lux = zeros(size(u,1),m);

                Qx = lx+A.'*v(:,i+1);
                Qu = lu+B.'*v(:,i+1);
                Qxx = lxx+A.'*S*A;
                Quu = luu+B.'*S*B;
                Qux = lux+B.'*S*A;

                Quur = Quu;

                d(:,i) = -Quur^(-1)*Qu;
                K(:,1+size(ds,1)*(i-1):size(ds,1)*i) = -Quur^(-1)*Qux;
    
                % Backward pass recursion
                S = Qxx+K(:,1+size(ds,1)*(i-1):size(ds,1)*i).'*(Quu*K(:,1+size(ds,1)*(i-1):size(ds,1)*i)+Qux)+Qux.'*K(:,1+size(ds,1)*(i-1):size(ds,1)*i);

                v(:,i) = Qx+K(:,1+size(ds,1)*(i-1):size(ds,1)*i).'*(Quu*d(:,i)+Qu)+Qux.'*d(:,i);
                % S = A.'*S*(A-B*K(:,1+size(ds,1)*(i-1):size(ds,1)*i))+Q; 
                % v(:,i) = (A-B*K(:,1+size(ds,1)*(i-1):size(ds,1)*i)).'*v(:,i+1)-K(:,1+size(ds,1)*(i-1):size(ds,1)*i).'*R*u(:,i)+Q*Sc(i,n+1:n+m).';
            end
    
            for i = 1:size(Sc,1)-1
                % Linearization 
                if (i > 1)
                    Phi0 = reshape(Sc(i-1,n+m+1:end), [m m]);     % STM from 0 to k-1
                else
                    Phi0 = eye(m);                                % STM from 0 to k-1
                end
                Phi = reshape(Sc(i,n+m+1:end), [m m]);            % STM from 0 to k
                A(1:m,1:m) = Phi*Phi0^(-1);                                             % Relative dynamics linear model
    
                switch (mode)
                    case 'Discrete'
                        B = A*b;                                                        % Control input matrix
                    case 'Continuous' 
                        B = (dt/2)*(eye(m)+A)*b;                                        % Control input matrix
                end    
    
                % Backward pass
                % du(:,i) = -(K(:,1+size(ds,1)*(i-1):size(ds,1)*i)*ds(:,i)+Ku(:,1+size(u,1)*(i-1):size(u,1)*i)*u(:,i)+Kv(:,1+size(v,1)*(i-1):size(v,1)*i)*v(:,i+1));
                du(:,i) = K(:,1+size(ds,1)*(i-1):size(ds,1)*i)*ds(:,i)+d(:,i);
                ds(:,i+1) = A*ds(:,i)+B*du(:,i);
            end
    
            % Forward pass        
            u = u+du;

            % APF
            
            switch (mode)
                case 'Discrete'
                    for i = 1:size(Sc,1)-1
                        U = zeros(1,n+m+m^2); 
                        U(n+4:n+6) = u(:,i).';
                        [~,s] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), [0 dt], Sc(i,:)+U, options);
                        Sc(i+1,:) = s(end,:);
                    end
    
                case 'Continuous'
                    dynamics = @(t,s)(vnlr_model(mu, u, t.', s.').');
                    [~, Sc, state] = MCPI([tspan(1) tspan(end)], tau, Sc, dynamics, N, 1e-12);
            end
    
            % Cost function evaluation
            Vp = 0; 
            switch (mode)
                case 'Discrete'
                    for i = 1:length(tspan)
                        Vp = Vp+0.5*Sc(i,n+1:n+m)*Q*Sc(i,n+1:n+m).'+0.5*u(:,i).'*R*u(:,i);
                    end
                case 'Continuous'           
                    for i = 1:length(tspan)
                        Vp = Vp+0.5*Sc(i,n+1:n+m)*Q*Sc(i,n+1:n+m).'+0.5*dt*u(:,i).'*R*u(:,i);
                    end
            end
    
            % Inequality penalty evaluations
            C = sqrt(dot(u,u,1))-0.01;
            Vp = Vp+dot(lambda,C)+0.5*rho*dot(C,C);

            % Convergence check 
            dV = abs(Vp-V);
            if (dV < tol(2)) 
                GoOn(2) = false;
            else
                iter(2) = iter(2)+1; 
                V = Vp;
            end
        end

        % Augmented Langrangian update 
        GoOn(1) = false; 
        max(C)
        if (max(C) < tol(1))
            GoOn(1) = false; 
        else
            iter(1) = iter(1)+1;
            lambda = max(0, lambda + rho .* C);
            rho = 2*rho;
        end
    end

    % Final output
    state.State = ~GoOn(1) && ~GoOn(2);
    state.Iterations = iter; 
    state.Error = dV;
end