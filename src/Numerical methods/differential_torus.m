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

% Outputs: - vector xf, converged/last corrected trajectory as well as other related parameter such as the TOF.
%          - boolean state, outputing the convergence of the algorithm.

% Methods: 

% New versions: select how many passes are allowed before stopping the
%               integration (the p-the Poincare image)

function [xf, state] = differential_torus(algorithm, mu, seed, maxIter, tol, nodes, constraint)
    %Implement the selected scheme 
    switch (algorithm)
        case 'Single shooting energy'
            [xf, state] = ssenergy_scheme(mu, seed, maxIter, tol, nodes, constraint);
        otherwise
            error('No valid option was selected');
    end
end

%% Auxiliary functions 
%Single shooting algorithm for constrained Jacobi constant 
function [xf, state]= ssenergy_scheme(mu, seed, maxIter, tol, nodes, constraint)
    %Constants 
    m = 6;                      %Phase space dimension 
    epsilon = 1e-3;             %Initial perturbation
    period = seed.Period;       %Periodic orbit period
    seed = seed.Trajectory;     %Periodic orbit
    
    %Preallocation 
    X = zeros(m*nodes+2,maxIter);                               %Preallocation of the free variables vector
    
    %Initialization of the free variables vector
    theta = 1/nodes*(0:2*pi:2*pi*(nodes-1));                    %Parametrization of the invariant curve
    refState = seed(end,:);                                     %Fixed point in the stroboscopic map
    Monodromy = reshape(refState(m+1:end), [m m]);              %Initial monodromy matrix
    [W, lambda] = eig(Monodromy);                               %Eigenspectrum of the monodromy matrix
    index = imag(lambda) ~= 0;                                  %Search for the center eigenvector
    lambda_c = lambda(index);                                   %Center eigenvalue
    v_c = W(:,any(index) == 1);                                 %Center eigenvectors
    X(end-1,1) = period;                                        %Initial orbital period
    X(end,1) = atan2(imag(lambda_c(1)), real(lambda_c(1)));     %Initial rotation angle
    
    %Compute the initial nodes 
    for i = 1:nodes
        X(1+m*(i-1):m*i,1) = refState(1:6).' + epsilon*(real(v_c(:,1))*cos(theta(i))-imag(v_c(:,1))*sin(theta(i)));
    end
    
    %Define the wavenumber on the domain 
    if (mod(nodes,2) == 0)
        k = -nodes/2:nodes/2;
    else
        k = -(nodes-1)/2:(nodes-1)/2;
    end
    
    %Define the Fourier transform operator 
    D = (1/nodes)*exp(-1i*k.'*theta);
    
    %Integration set up
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Preallocation of the error and the sensibility matrices 
    urot = zeros(nodes,m);                                      %First return of the nodes without the rotation
    e = zeros(m*nodes+1, 1);                                    %Error vector
    bSTM = zeros(m*nodes,m*nodes);                              %Augmented state transition matrix
    dJ = zeros(1,nodes*m);                                      %Jacobian matrix of the Jacobi constant
    dF = zeros(m*nodes, 1);                                     %Time derivative of the nodes
    
    %Set up differential correction scheme
    GoOn = true;                %Convergence flag
    iter = 1;                   %Initial iteration
    
    %Differential corrector
    while (GoOn) && (iter < maxIter)
        %Set up the integration 
        T = real(X(end-1,iter));        %Strobocopic time
        tspan = 0:dt:T;                 %Integration time
        
        %Reinitiate the Jacobi constant 
        J = 0; 
        
        %Compute the rotation operator
        Q = diag(exp(-1i*k*X(end,iter)));                   %Rotation eigenvalues
        dQ = diag(k)*diag(-1i*exp(-1i*k*X(end,iter)));      %Derivative of the application with respect to the rotation angle  
        R = D^(-1)*Q*D;                                     %Rotation operation
        
        %Perform the integration 
        for i = 1:nodes
            S0 = X(1+m*(i-1):m*i,iter);         %State vector 
            STM = reshape(eye(m), [m^2 1]);     %Initial state transition matrix
            s0 = [S0; STM];                     %Complete initial conditions
            
            %Integration
            [~, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, s0, options);  
            
            %Compute the associated Jacobi Constant 
            J = J + jacobi_constant(mu, S(end,1:6).');
            
            %Desrotate the state 
            urot(i,1:m) = S(end,1:6);
            
            %Time vector field and associated STM
            dF(1+m*(i-1):m*i,:) = cr3bp_equations(mu, direction, false, 0, S(end,1:6).');
            bSTM(1+m*(i-1):m*i,1+m*(i-1):m*i) = reshape(S(end,m+1:end), [m m]);
        end
        
        %Desrotate the states 
        u = R*urot;
        
        %Compute the error vector 
        for i = 1:nodes
            e(1+m*(i-1):m*i,1) = u(i,:).'-X(1+m*(i-1):m*i,iter);       %Rotation constraint
        end
        e(end,1) = (1/nodes)*J-constraint;                             %Energy constraint
        
        %Compute the sensibility matrix
        DG = kron(R,eye(m))*bSTM-eye(m*nodes);  %STM-like sensibility matrix
        dRho = D^(-1)*dQ*D*u;                   %Derivative with respect to the rotation angle
        dRho = reshape(dRho, [m*nodes 1]);      %Derivative with respect to the rotation angle
        
        for i = 1:nodes
            dJ(1,1+m*(i-1):m*i) = (1/nodes)*jacobi_gradient(mu, u(i,:).');
        end
        
        A = [DG dF dRho; dJ zeros(1,2)];
        
        %Compute the variation 
        ds = -real(A\e); 
        norm(ds)
                        
        %Convergence analysis 
        if (norm(ds) < tol)
            GoOn = false;               %Finish the correction process
        else
            X(:,iter+1) = X(:,iter)+ds; %Update the variables vector 
            iter = iter+1;              %Update the iterations
        end
    end

    %Output
    state.State = ~GoOn;                %Convergence flag
    state.Error = 1;                    %Finall correction error
    state.Iter = iter;                  %Number of needed iterations
    
end

