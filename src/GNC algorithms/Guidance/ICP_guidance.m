%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date:  09/07/22
% File: ICP_guidance.m 
% Issue: 0 
% Validated: 26/07/22

%% Impulsive Center Manifold Phasing %%
% This script contains the function to compute a phasing guidance trajectory using
% the center manifold expansion of the relative dynamics

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar tf, the final time of flight
%         - vector s0, initial conditions of both the target and the chaser
%         - scalar epsilon, the manifold displaement factor
%         - scalar tol, the differential corrector scheme tolerance for the
%           constrained maneuver

% Output: - array S, the rendezvous relative trajectory
%         - structure state, containing the results of the differential
%           corrector

% New versions: 

%% Auxiliary functions 
function [S, state, Tref] = ICP_guidance(mu, T, dtheta, K, s0, epsilon, tol)
    %Constants 
    m = 6;       % Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
    
    %Integration tolerance setup and time span
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); %Tolerances
    dt = 1e-3;                                             %Time step  
    tspan = 0:dt:T;                                        %Integration time span

    %Computation of the monodromy matrix
    S0 = [s0(1:m); s0(m+1:end)-s0(1:m)];                   %Initial relative chaser conditions
    Phi = eye(m);                                          %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                           %Reshape the initial STM
    S0 = [S0; Phi];                                        %Complete phase space + linear variational initial conditions

    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, S0, options); 
        
    %Dimensional analysis of the center manifold of interest
    Monodromy = reshape(Sn(end,2*m+1:end), [m m]);         %Initial monodromy matrix
    [~, lambda] = eig(Monodromy);                          %Eigenspectrum of the monodromy matrix
    lambda = diag(lambda);                                 %Consider only the pure eigenalues
    index = imag(lambda) ~= 0;                             %Search for the center eigenvector
    lambda_c = lambda(index);                              %Center eigenvalue
    dim = nnz(lambda_c);                                   %Dimension of the relative center manifold

    % Final algorithm
    if (dim == 2)
        [S, state, Tref] = one_torus(mu, T, dtheta, K, s0, epsilon, tol);
    elseif (dim == 4)
        [S, state, Tref] = two_torus(mu, T, dtheta, K, s0, epsilon, tol);
    else
        warning('No relative center manifold exists');
        S = []; 
        state.State = false; 
        Tref = T; 
    end
end

%% Auxiliary functions 
% Phasing on a 2D center manifold 
function [S, state, Tref] = one_torus(mu, T, dtheta, K, s0, epsilon, tol)
    %Constants 
    m = 6;       % Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end

    %Compute the reference period 
    Tref = T*(1-dtheta/(2*pi*K));                               %Needed relative period
    
    %Integration time span
    dt = 1e-3;                                                  %Time step  
    tspan = 0:dt:T;                                             %Integration time span

    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 

    %Computation of the monodromy matrix
    s0 = [s0(1:m); s0(m+1:end)-s0(1:m)];                   %Initial relative chaser conditions
    Phi = eye(m);                                          %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                           %Reshape the initial STM
    s0 = [s0; Phi];                                        %Complete phase space + linear variational initial conditions

    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 
        
    %Differential corrector setup
    GoOn = true;                                           %Convergence boolean
    maxIter = 100;                                         %Maximum number of iterations
    iter = 1;                                              %Initial iteration

    %Initial constants and parameters
    nodes = 2*K+1;                                         %Number of points on the invariant curve
    theta = 2*pi*(0:(nodes-1))/nodes;                      %Parametrization of the invariant curve
    Monodromy = reshape(Sn(end,2*m+1:end), [m m]);         %Initial monodromy matrix
    [W, lambda] = eig(Monodromy);                          %Eigenspectrum of the monodromy matrix
    lambda = diag(lambda);                                 %Consider only the pure eigenalues
    index = imag(lambda) ~= 0;                             %Search for the center eigenvector
    lambda_c = lambda(index);                              %Center eigenvalue
    v_c = W(:,index);                                      %Center eigenvector
    X = zeros(m*nodes+2,maxIter);                          %Preallocation of the free variables state vector

    %Definition of the phasing period and rotation number
    X(end-1,1) = T;                                             %Initial orbital period or stroboscopic time
    X(end,1) = atan2(imag(lambda_c(1)), real(lambda_c(1)));     %Initial rotation angle

    %Initial quasi-periodic orbit conditions
    for i = 1:nodes
        X(1+m*(i-1):m*i,1) = epsilon*(real(v_c(:,1))*cos(theta(i))-imag(v_c(:,1))*sin(theta(i)));
    end

    %Define the wavenumber on the domain 
    if (mod(nodes,2) == 0)
        k = -nodes/2:nodes/2;
    else
        k = -(nodes-1)/2:(nodes-1)/2;
    end
    
    %Define the Fourier transform operator and its inverse kernel
    D = exp(-1i*k.'*theta)/nodes;
    invD = D^(-1);

    %Preallocation of the error and the sensibility matrices 
    urot = zeros(nodes,m);                              %First return of the nodes without the rotation
    bSTM = zeros(m*nodes,m*nodes);                      %Augmented state transition matrix
    dF = zeros(m*nodes, 1);                             %Time derivative of the nodes

    %Minimum velocity impulse 
    Vmin = 1e-7;
        
    while ((GoOn) && (iter < maxIter))
        %Torus variables
        T = X(end-1,iter);                              %Strobocopic time
        rho = X(end,iter);                              %Rotation operator
        tspan = 0:dt:T;                                 %Integration time

        Q = diag(exp(-1i*k*rho));                       %Rotation eigenvalues
        dQ = diag(k)*diag(-1i*exp(-1i*k*rho));          %Derivative of the application with respect to the rotation angle  
        R = invD*Q*D;                                   %Rotation operation

        %Perform the stroboscopic mapping
        for i = 1:nodes
            %Complete relative initial conditions
            S0 = [s0(1:m)-s0(m+1:2*m); X(1+m*(i-1):m*i,iter); s0(2*m+1:end)];            
            
            %Integration
            [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, S0, options);  
                        
            %First return to the invariant curve states
            urot(i,1:m) = Saux(end,m+1:2*m);
            
            %Time vector field and associated STM
            f = nlr_model(mu, true, false, false, 'Encke', T, Saux(end,1:2*m).');
            dF(1+m*(i-1):m*i,:) = f(m+1:2*m);
            bSTM(1+m*(i-1):m*i,1+m*(i-1):m*i) = reshape(Saux(end,2*m+1:end), [m m]);
        end

        %Minimum impulse constraint 
        dv = zeros(nodes,1);
        VA = zeros(nodes,m*nodes);
        for i = 1:nodes
            index = 4+m*(i-1);
            dv(i) = norm(X(index:index+2,iter))-Vmin;
            VA(i,1+m*(i-1):m*i) = [zeros(1,3) X(index:index+2,iter).'/norm(X(index:index+2,iter))];
        end

        %Compute the error vector 
        u = R*urot;
        error = [reshape(u.', [nodes*m,1])-X(1:nodes*m,iter); T-Tref; dv];

        %Compute the sensibility matrix
        DG = kron(R,eye(m))*bSTM-eye(m*nodes);  %STM-like sensibility matrix
        dRho = invD*dQ*D*urot;                  %Derivative with respect to the rotation angle
        dRho = reshape(dRho, [m*nodes 1]);      %Derivative with respect to the rotation angle

        %Compute the Newton-Rhapson update
        A = [DG dF dRho; zeros(1,m*nodes) 1 0; VA zeros(nodes,2)]; %Complete sensibility matrix
        ds = real(-pinv(A)*error);                                 %Newton-Rhapson innovation
        X(:,iter+1) = X(:,iter)+ds;                                %Update the variables vector 

        %Convergence analysis 
        if (norm(ds) < tol)
            GoOn = false;                       %Finish the correction process
        else
            iter = iter+1;                      %Update the iterations
        end
    end
     
    %Final invariant curve initial conditions
    ds = reshape(X(1:end-2,iter), [nodes m]);               %Final state variation

    %Select the final closest one in position space
    dS = sqrt(dot(ds(:,1:3),ds(:,1:3),2));
    [~, index] = sort(dS, 'descend'); 
    ds = ds(index,:);
    ds = ds(end-K:end,:);

    %Final mappings 
    S.Trajectory = [repmat(s0(1:m).', [K+1 1]) repmat(s0(m+1:2*m).', [K+1 1])+ds];               

    S.Period = X(end-1,iter);                               %Final stroboscopic time
    S.Rotation = X(end,iter);                               %Final rotation number
    state.State = ~GoOn;                                    %Final convergence flag 
    state.Iter = iter;                                      %Final iteration 
    state.Error = norm(error);                              %Final error
end

% Phasing on a 4D center manifold 
function [S, state] = two_torus(mu, T, dtheta, K, s0, epsilon, tol)
    %Constants 
    m = 6;       % Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end

    %Compute the reference period 
    Tref = T*(1-dtheta/(2*pi*K));                               %Needed relative period
    
    %Integration time span
    dt = 1e-3;                                                  %Time step  
    tspan = 0:dt:T;                                             %Integration time span

    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 

    %Computation of the monodromy matrix
    s0 = [s0(1:m); s0(m+1:end)-s0(1:m)];                   %Initial relative chaser conditions
    Phi = eye(m);                                          %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                           %Reshape the initial STM
    s0 = [s0; Phi];                                        %Complete phase space + linear variational initial conditions

    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 
        
    %Differential corrector setup
    GoOn = true;                                           %Convergence boolean
    maxIter = 50;                                          %Maximum number of iterations
    iter = 1;                                              %Initial iteration

    %Initial constants and parameters
    nodes = 2*K+1;                                         %Number of points on the invariant curve
    theta = 2*pi*(0:(nodes-1))/nodes;                      %Parametrization of the invariant curve
    Monodromy = reshape(Sn(end,2*m+1:end), [m m]);         %Initial monodromy matrix
    [W, lambda] = eig(Monodromy);                          %Eigenspectrum of the monodromy matrix
    lambda = diag(lambda);                                 %Consider only the pure eigenalues
    index = imag(lambda) ~= 0;                             %Search for the center eigenvector
    lambda_c = lambda(index);                              %Center eigenvalue
    v_c = W(:,index);                                      %Center eigenvector
    X = zeros(m*nodes+2,maxIter);                          %Preallocation of the free variables state vector

    %Definition of the phasing period and rotation number
    X(end-1,1) = T;                                             %Initial orbital period or stroboscopic time
    X(end,1) = atan2(imag(lambda_c(1)), real(lambda_c(1)));     %Initial rotation angle

    %Initial quasi-periodic orbit conditions
    epsilon = 1e-4;                                             % Center manifold displacement
    for i = 1:nodes
        X(1+m*(i-1):m*i,1) = epsilon*(real(v_c(:,1))*cos(theta(i))-imag(v_c(:,1))*sin(theta(i)));
    end

    %Define the wavenumber on the domain 
    if (mod(nodes,2) == 0)
        k = -nodes/2:nodes/2;
    else
        k = -(nodes-1)/2:(nodes-1)/2;
    end
    
    %Define the Fourier transform operator and its inverse kernel
    D = exp(-1i*k.'*theta)/nodes;
    invD = D^(-1);

    %Preallocation of the error and the sensibility matrices 
    urot = zeros(nodes,m);                              %First return of the nodes without the rotation
    bSTM = zeros(m*nodes,m*nodes);                      %Augmented state transition matrix
    dF = zeros(m*nodes, 1);                             %Time derivative of the nodes
        
    while ((GoOn) && (iter < maxIter))
        %Torus variables
        T = X(end-1,iter);                              %Strobocopic time
        rho = X(end,iter);                              %Rotation operator
        tspan = 0:dt:T;                                 %Integration time

        Q = diag(exp(-1i*k*rho));                       %Rotation eigenvalues
        dQ = diag(k)*diag(-1i*exp(-1i*k*rho));          %Derivative of the application with respect to the rotation angle  
        R = invD*Q*D;                                   %Rotation operation

        %Perform the stroboscopic mapping
        for i = 1:nodes
            %Complete relative initial conditions
            S0 = [s0(1:m)-s0(m+1:2*m); X(1+m*(i-1):m*i,iter); s0(2*m+1:end)];            
            
            %Integration
            [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, S0, options);  
                        
            %First return to the invariant curve states
            urot(i,1:m) = Saux(end,m+1:2*m);
            
            %Time vector field and associated STM
            f = nlr_model(mu, true, false, false, 'Encke', T, Saux(end,1:2*m).');
            dF(1+m*(i-1):m*i,:) = f(m+1:2*m);
            bSTM(1+m*(i-1):m*i,1+m*(i-1):m*i) = reshape(Saux(end,2*m+1:end), [m m]);
        end

        %Compute the error vector 
        u = R*urot;
        error = [reshape(u.', [nodes*m,1])-X(1:nodes*m,iter); T-Tref];

        %Compute the sensibility matrix
        DG = kron(R,eye(m))*bSTM-eye(m*nodes);  %STM-like sensibility matrix
        dRho = invD*dQ*D*urot;                  %Derivative with respect to the rotation angle
        dRho = reshape(dRho, [m*nodes 1]);      %Derivative with respect to the rotation angle

        %Compute the Newton-Rhapson update
        A = [DG dF dRho; zeros(1,m*nodes) 1 0]; %Complete sensibility matrix
        ds = real(-pinv(A)*error);              %Newton-Rhapson innovation
        X(:,iter+1) = X(:,iter)+ds;             %Update the variables vector 

        %Convergence analysis 
        if (norm(ds) < tol)
            GoOn = false;                       %Finish the correction process
        else
            iter = iter+1;                      %Update the iterations
        end
    end
     
    %Final invariant curve initial conditions
    ds = reshape(X(1:end-2,iter), [nodes m]);               %Final state variation

    %Select the final closest one in position space
    dS = sqrt(dot(ds(:,1:3),ds(:,1:3),2));
    [~, index] = sort(dS, 'descend'); 
    ds = ds(index,:);
    ds = ds(end-K:end,:);

    %Final mappings 
    S.Trajectory = [repmat(s0(1:m).', [K+1 1]) repmat(s0(m+1:2*m).', [K+1 1])+ds];               

    S.Period = X(end-1,iter);                               %Final stroboscopic time
    S.Rotation = X(end,iter);                               %Final rotation number
    state.State = ~GoOn;                                    %Final convergence flag 
    state.Iter = iter;                                      %Final iteration 
    state.Error = norm(error);                              %Final error
end