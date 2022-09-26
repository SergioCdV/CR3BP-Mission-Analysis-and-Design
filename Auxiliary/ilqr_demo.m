%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 24/09/22 % 

%% GNC 3: iLQR control law %% 
% This script provides an interface to test iLQR rendezvous strategies for
% rendezvous missions

% Units are non-dimensional and solutions are expressed in the synodic
% reference frame as defined by Howell, 1984

%% Set up %%
% Set up graphics 
set_graphics();

% Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
% Phase space dimension 
n = 6; 

% Time span 
dt = 1e-3;                                      % Time step
tf = 10;                                         % Rendezvous time

%% GNC algorithms definition 
GNC.Control.LQR.Q = 1e0*eye(2);                 % Penalty on the state error
GNC.Control.LQR.M = 0.8*1e0*eye(1);                 % Penalty on the control effort

% iLQR control law 
int = 0;                                        % Integral of the relative position
slqr0 = [10 1];                                 % Initial conditions

% Compute the trajectory
[tspan, St, u, state] = iLQR_control_OS(tf, slqr0, GNC); 

%% Results %% 
figure
plot(tspan,St);


%Plot results 
% figure
% view(3) 
% hold on
% plot3(Sn(:,1), Sn(:,2), Sn(:,3)); 
% plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
% hold off
% legend('Target motion', 'Chaser motion'); 
% xlabel('Synodic $x$ coordinate');
% ylabel('Synodic $y$ coordinate');
% zlabel('Synodic $z$ coordinate');
% grid on;
% title('Reconstruction of the natural chaser motion');
% 
% %Plot relative phase trajectory
% figure
% view(3) 
% plot3(St(:,7), St(:,8), St(:,9)); 
% xlabel('Synodic $x$ coordinate');
% ylabel('Synodic $y$ coordinate');
% zlabel('Synodic $z$ coordinate');
% grid on;
% title('Motion in the relative configuration space');

%% Auxiliary functions
function [tspan, Sc, u, state] = iLQR_control_DI(tf, s0, GNC)
    % Constants of the model 
    tol = 1e-5;                   % Convergence tolerance

    B = [0; 1];                   % Control input matrix
    A = [0 1; 0 0];               % PID-structured Jacobian matrix

    % Controller definition
    Q = GNC.Control.LQR.Q;                          % Penalty matrix on the state error
    R = GNC.Control.LQR.M;                          % Penalty matrix on the control effort

    % Generate the initial LQR/SDRE guess 
    dt = 1e-3;                                          % Nondimensional time step
    tspan = 0:dt:tf;                                    % Integration time span

    % Preallocation 
    Sc = repmat(s0,length(tspan),1);                    % Nominal trajectory
    u = zeros(1,length(tspan));                         % Nominal control law
    du = zeros(1,size(Sc,1));                           % Update to the control law 
    ds = zeros(length(s0),size(Sc,1));                  % Update to the rendezvous trajectory
    v = zeros(length(s0),size(Sc,1));                   % Feedforward term to the rendezvous trajectory

    K = zeros(size(u,1),size(ds,1)*size(Sc,1));         % LQR gain
    Ku = zeros(size(u,1),size(u,1)*size(Sc,1));         % Feedback gain
    Kv = zeros(size(u,1),size(v,1)*size(Sc,1));         % Feeforward gain
    
    % iLQR iterative loop 
    GoOn = true;                            % Convergence boolean 
    maxIter = 100;                          % Maximum number of iterations
    iter = 1;                               % Initial iteration 

    for i = 1:size(Sc,1)-1
        F(:,1) = [Sc(i,2); u(:,i)];
        aux = Sc(i,:).'+F(:,1)*dt;
        F(:,2) = [aux(2); u(:,i+1)];
        Sc(i+1,:) = Sc(i,:) + dt*sum(F,2).'/2;
    end

    % Boundary conditions 
    ds(:,1) = zeros(length(s0),1);              % Initial update

    V = 0; 
    for i = 1:length(tspan)
        V = V+0.5*Sc(i,:)*Q*Sc(i,:).'+0.5*u(:,i).'*R*u(:,i);
    end

    while (GoOn && iter < maxIter)
        S = Q;                                  % Penalty weight matrix final conditions
        v(:,end) = Q*Sc(end,:).';               % Final feedforward term
    
        for i = size(Sc,1)-1:-1:1
            % Controller gains
            K(:,1+size(ds,1)*(i-1):size(ds,1)*i) = (B.'*S*B+R)^(-1)*B.'*S*A;
            Kv(:,1+size(v,1)*(i-1):size(v,1)*i) = (B.'*S*B+R)^(-1)*B.';
            Ku(:,1+size(u,1)*(i-1):size(u,1)*i) = (B.'*S*B+R)^(-1)*R;

            % Backward pass recursion
            S = A.'*S*(A-B*K(:,1+size(ds,1)*(i-1):size(ds,1)*i))+Q; 
            v(:,i) = (A-B*K(:,1+size(ds,1)*(i-1):size(ds,1)*i)).'*v(:,i+1)-K(:,1+size(ds,1)*(i-1):size(ds,1)*i).'*R*u(:,i)+Q*Sc(i,:).';
        end

        for i = 1:size(Sc,1)-1
            % Backward pass
            du(:,i) = -(K(:,1+size(ds,1)*(i-1):size(ds,1)*i)*ds(:,i)+Ku(:,1+size(u,1)*(i-1):size(u,1)*i)*u(:,i)+Kv(:,1+size(v,1)*(i-1):size(v,1)*i)*v(:,i+1));

            % Forward pass        
            u(:,i) = u(:,i)+du(:,i);
            ds(:,i+1) = A*ds(:,i)+B*du(:,i);
        end

        % Rollout
        for i = 1:size(Sc,1)-1
            F(:,1) = [Sc(i,2); u(:,i)];
            aux = Sc(i,:).'+F(:,1)*dt;
            F(:,2) = [aux(2); u(:,i+1)];
            Sc(i+1,:) = Sc(i,:) + dt*sum(F,2).'/2;
        end
        
        Vp = 0; 
        for i = 1:length(tspan)
            Vp = Vp+0.5*Sc(i,:)*Q*Sc(i,:).'+0.5*u(:,i).'*R*u(:,i);
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

function [tspan, Sc, u, state] = iLQR_control_OS(tf, s0, GNC)
    % Constants of the model 
    tol = 1e-5;                   % Convergence tolerance

    B = [0; 1];                   % Control input matrix
    A = [0 1; -1 0];              % PID-structured Jacobian matrix

    % Controller definition
    Q = GNC.Control.LQR.Q;                          % Penalty matrix on the state error
    R = GNC.Control.LQR.M;                          % Penalty matrix on the control effort

    % Generate the initial LQR/SDRE guess 
    dt = 1e-3;                                          % Nondimensional time step
    tspan = 0:dt:tf;                                    % Integration time span

    % Preallocation 
    Sc = repmat(s0,length(tspan),1);                    % Nominal trajectory
    u = zeros(1,length(tspan));                         % Nominal control law
    du = zeros(1,size(Sc,1));                           % Update to the control law 
    ds = zeros(length(s0),size(Sc,1));                  % Update to the rendezvous trajectory
    v = zeros(length(s0),size(Sc,1));                   % Feedforward term to the rendezvous trajectory

    K = zeros(size(u,1),size(ds,1)*size(Sc,1));         % LQR gain
    Ku = zeros(size(u,1),size(u,1)*size(Sc,1));         % Feedback gain
    Kv = zeros(size(u,1),size(v,1)*size(Sc,1));         % Feeforward gain
    
    % iLQR iterative loop 
    GoOn = true;                            % Convergence boolean 
    maxIter = 100;                          % Maximum number of iterations
    iter = 1;                               % Initial iteration 

    for i = 1:size(Sc,1)-1
        F(:,1) = [Sc(i,2); -sin(Sc(i,1))^2+u(:,i)];
        aux = Sc(i,:).'+F(:,1)*dt;
        F(:,2) = [aux(2); -sin(aux(1))^2+u(:,i+1)];
        Sc(i+1,:) = Sc(i,:) + dt*sum(F,2).'/2;
    end

    % Boundary conditions 
    ds(:,1) = zeros(length(s0),1);              % Initial update

    V = 0; 
    for i = 1:length(tspan)
        V = V+0.5*Sc(i,:)*Q*Sc(i,:).'+0.5*u(:,i).'*R*u(:,i);
    end

    while (GoOn && iter < maxIter)
        S = Q;                                  % Penalty weight matrix final conditions
        v(:,end) = Q*Sc(end,:).';               % Final feedforward term
   
        for i = size(Sc,1)-1:-1:1
            A(2,1) = -2*sin(Sc(i,1))*cos(Sc(i,1));

            % Controller gains
            K(:,1+size(ds,1)*(i-1):size(ds,1)*i) = (B.'*S*B+R)^(-1)*B.'*S*A;
            Kv(:,1+size(v,1)*(i-1):size(v,1)*i) = (B.'*S*B+R)^(-1)*B.';
            Ku(:,1+size(u,1)*(i-1):size(u,1)*i) = (B.'*S*B+R)^(-1)*R;

            % Backward pass recursion
            S = A.'*S*(A-B*K(:,1+size(ds,1)*(i-1):size(ds,1)*i))+Q; 
            v(:,i) = (A-B*K(:,1+size(ds,1)*(i-1):size(ds,1)*i)).'*v(:,i+1)-K(:,1+size(ds,1)*(i-1):size(ds,1)*i).'*R*u(:,i)+Q*Sc(i,:).';
        end

        for i = 1:size(Sc,1)-1
            A(2,1) = -2*sin(Sc(i,1))*cos(Sc(i,1));

            % Backward pass
            du(:,i) = -(K(:,1+size(ds,1)*(i-1):size(ds,1)*i)*ds(:,i)+Ku(:,1+size(u,1)*(i-1):size(u,1)*i)*u(:,i)+Kv(:,1+size(v,1)*(i-1):size(v,1)*i)*v(:,i+1));

            % Forward pass        
            u(:,i) = u(:,i)+du(:,i);
            ds(:,i+1) = A*ds(:,i)+B*du(:,i);
        end

        % Rollout
        for i = 1:size(Sc,1)-1
            F(:,1) = [Sc(i,2); -sin(Sc(i,1))^2+u(:,i)];
            aux = Sc(i,:).'+F(:,1)*dt;
            F(:,2) = [aux(2); -sin(aux(1))^2+u(:,i+1)];
            Sc(i+1,:) = Sc(i,:) + dt*sum(F,2).'/2;
        end
        
        Vp = 0; 
        for i = 1:length(tspan)
            Vp = Vp+0.5*Sc(i,:)*Q*Sc(i,:).'+0.5*u(:,i).'*R*u(:,i);
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