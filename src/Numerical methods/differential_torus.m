%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 03/01/21
% File: differential_correction.m 
% Issue: 0 
% Validated: 

%% Transfer differential correction %%
% This function contains the algorithm to compute several correction
% algorithms to generate low-budget transfers in the C3RBP context.

% Inputs: - string algorithm, selecting the differential correction scheme to use. 
%         - double mu, the reduced gravitational parameter of the system.
%         - object parking_orbit, defining the initial parking orbit.
%         - object target_orbit, which vary depending on the algorithm in use. 
%         - int maxIter, number of maximum allowed corrections.
%         - double tol, tolerance to stop the correction. 

% Outputs: - structure S, converged/last corrected initial conditions as well as other related parameter such as the TOF.
%          - boolean state, outputing the convergence of the algorithm.

% Methods: 

% New versions: select how many passes are allowed before stopping the
%               integration (the p-the Poincare image)

function [S, state] = differential_torus(algorithm, mu, seed, tol)
    %Implement the selected scheme 
    switch (algorithm)
        case 'Single shooting energy'
            L = libration_points(mu);
            L0 = [L(1:3,1).' zeros(1,3)];
            [S, state] = differential_rtorus(mu, seed.Period, [L0 seed.Trajectory(1,1:6)-L0], tol);
        otherwise
            error('No valid option was selected');
    end
end

%% Auxiliary functions 
%Single shooting algorithm for constrained Jacobi constant 
function [S, state]= ssenergy_scheme(mu, seed, maxIter, tol, nodes)
    %Constants 
    m = 6;                      %Phase space dimension 
    period = seed.Period;       %Periodic orbit period
    seed = seed.Trajectory;     %Periodic orbit

    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
        
    %Initial constants and parameters
    theta = 2*pi*(0:(nodes-1))/nodes;                      %Parametrization of the invariant curve
    Monodromy = reshape(seed(end,m+1:end), [m m]);         %Initial monodromy matrix
    [W, lambda] = eig(Monodromy);                          %Eigenspectrum of the monodromy matrix
    lambda = diag(lambda);                                 %Consider only the pure eigenalues
    index = imag(lambda) ~= 0;                             %Search for the center eigenvector
    lambda_c = lambda(index);                              %Center eigenvalue
    v_c = W(:,index);                                      %Center eigenvector
    X = zeros(m*nodes+2,maxIter);                          %Preallocation of the free variables state vector

    %Definition of the phasing period and rotation number
    X(end-1,1) = period;                                        %Initial orbital period or stroboscopic time
    X(end,1) = atan2(imag(lambda_c(1)), real(lambda_c(1)));     %Initial rotation angle
    
    %Initial quasi-periodic orbit conditions
    epsilon = 1e-3;                                             % Center manifold displacement
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
    
    %Integration set up                  
    dt = 1e-3;                                                  %Time step
       
    %Set up differential correction scheme
    GoOn = true;                %Convergence flag
    iter = 1;                   %Initial iteration
    
    %Preallocation of the error and the sensibility matrices 
    urot = zeros(nodes,m);                              %First return of the nodes without the rotation
    bSTM = zeros(m*nodes,m*nodes);                      %Augmented state transition matrix
    dF = zeros(m*nodes, 1);                             %Time derivative of the nodes
    dJ = zeros(1,nodes*m);                              %Jacobian matrix of the Jacobi constant
        
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
            S0 = [X(1+m*(i-1):m*i,iter); reshape(eye(m), [m^2 1])];            
            
            %Integration
            [~, Saux] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, S0, options);  
                        
            %First return to the invariant curve states
            urot(i,1:m) = Saux(end,m+1:2*m);
            
            %Time vector field and associated STM
            bSTM(1+m*(i-1):m*i,1+m*(i-1):m*i) = reshape(Saux(end,2*m+1:end), [m m]);
        end

        %Compute the error vector 
        u = R*urot;
        error = reshape(u.', [nodes*m,1])-X(1:nodes*m,iter);

        %Compute the sensibility matrix
        DG = kron(R,eye(m))*bSTM-eye(m*nodes);  %STM-like sensibility matrix
        dRho = invD*dQ*D*urot;                  %Derivative with respect to the rotation angle
        dRho = reshape(dRho, [m*nodes 1]);      %Derivative with respect to the rotation angle

        for i = 1:nodes 
            f = cr3bp_equations(mu, true, false, period, (seed(1,1:m)+u(i,:)).');
            dF(1+m*(i-1):m*i,:) = f(m+1:2*m);
        end

        %Compute the Newton-Rhapson update
        A = [DG dF dRho];                       %Complete sensibility matrix
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
    S.Trajectory = repmat(seed(1:m), nodes, 1)+reshape(X(1:end-2,iter), [nodes m]);     

    S.Period = X(end-1,iter);                               %Final stroboscopic time
    S.Rotation = X(end,iter);                               %Final rotation number
    state.State = ~GoOn;                                    %Final convergence flag 
    state.Iter = iter;                                      %Final iteration 
    state.Error = norm(error);                              %Final error
end

