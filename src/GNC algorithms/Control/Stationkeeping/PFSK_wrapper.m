%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/11/21
% File: PFSK_wrapper.m 
% Issue: 0 
% Validated: 30/11/21

%% Primer Floquet Stationkeeping %%
% This script contains the function to compute the control law by means of an PFSK controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar T, the target orbit reference period
%         - vector s0, initial conditions of the target spacecraft
%         - structure constraint, specifying any constraint on the maneuver
%         - string cost_function, to optimize on the l1 or l2 norm of the
%           control vector
%         - vector Tmax, the maximum available thrust

% Output: - array Sc, the stationkeeping trajectory
%         - array u, containing the required control law
%         - structure state, with the corresponding differential corrector
%           figures of merit

% New versions: 

function [S, u] = PFSK_wrapper(mu, T, tf, s0, constraint, cost_function, Tmax)
    %Constants 
    m = 6;       %Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
    
    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    dt = 1e-3;                                                  %Integration time step  
    
    %Initial conditions and integration
    Phi = eye(m);                                               %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                                %Reshape the initial STM
    s0 = [s0; Phi];                                             %Complete phase space + linear variational initial conditions
    
    tspan = 0:dt:tf;                                            %Integration time span
    [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);

    %Energy constraint
    constraint_flag = constraint.Flag;                          %Constraint flag
    if (constraint_flag)
        Cref = constraint.JacobiReference;                      %Reference Jacobi Constant
    end

    %Definition of the GNC structure
    GNC.Algorithms.Guidance = '';                       %Guidance algorithm
    GNC.Algorithms.Navigation = '';                     %Navigation algorithm
    GNC.Algorithms.Control = 'PFSK';                    %Control algorithm
    
    GNC.Guidance.Dimension = 6;                         %Dimension of the guidance law
    GNC.Control.Dimension = 3;                          %Dimension of the control law
    
    GNC.System.mu = mu;                                 %Systems's reduced gravitational parameter
    GNC.Control.PFSK.CostFunction = cost_function;      %Cost function to optimize
    GNC.Control.PFSK.MaxThrust = Tmax;                  %Maximum available thurst

    %Floquet analysis                                             
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), 0:dt:T, s0, options);

    Monodromy = reshape(Sn(end,2*m+1:end), [m,m]);              %Monodromy matrix of the relative orbit
    [F,J] = eig(Monodromy);                                     %Eigenspectrum of the monodromy matrix
    J = diag(log(diag(J))/T);                                   %Floquet exponents
    for j = 1:size(F,2)
        F(:,j) = F(:,j)/J(j,j);                                 %Floquet basis initial conditions
    end

    GNC.Control.PFSK.FloquetExponents = J;                      %Floquet exponents of the reference trajectory
    GNC.Control.PFSK.FloquetDirections = F;                     %Floquet directions of the reference trajectory

    %Compute the initial guess
    M = reshape(Saux(end,2*m+1:end), [m m]);                    %Instantenous Monodromy matrix
    P = F*expm(-J*mod(tf,T));                                   %Full Floquet projection matrix
    E = M*P;                                                    %Instantenous Floquet projection matrix

    alpha(:,1) = F^(-1)*Saux(1,m+1:2*m).';                      %Initial Floquet coordinates
    alpha(:,2) = E^(-1)*Saux(end,m+1:2*m).';                    %Final Floquet coordinates
    lambda(:,2) = [1; -1; zeros(4,1)];                          %Final co-state guess 
    lambda(:,1) = expm(-J.'*tf)*lambda(:,2);                    %Initial co-state guess
    GNC.Control.PFSK.InitialPrimer = lambda(:,1);               %Initial primer vector

    %Initial guess for the numerical method
    nsteps = length(0:dt:tf);
    solinit.x = linspace(0, tf, nsteps);
    solinit.y = repmat([alpha(:,1); lambda(:,1); s0], 1, nsteps);

    %Set optimizer options
    tol = 1e-10;
    options = bvpset('RelTol', tol, 'AbsTol', tol*ones(60,1), 'Nmax', 2000);

    %Solve the problem 
    sol = bvp4c(@(t,s)dynamics(mu, t, s, T, GNC), @(x,x2)bounds(x, x2, alpha(:,1)), solinit, options);

    %Output 
    Saux = [sol.x sol.y];                                   %Final trajectory in the Floquet space
    GNC.Control.PFSK.InitialPrimer = Saux(1,m+1:2*m);       %Initial primer vector

    %Re-integration of the trajectory
    [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan, s0, options);

    %Re-assembly of the control vector 
    [~, ~, u] = GNCc_handler(GNC, S(:,1:m), S(:,m+1:end), tspan);
end

%% Auxiliary functions 
% Optimal problem dynamics 
function [ds] = dynamics(mu, t, s, T, GNC)
    % Constants 
    m = 6;          %Phase space dimension 

    % Physical phase space dynamics 
    [ds] = nlr_model(mu, true, false, true, 'Encke', t, s(2*m+1:end), GNC);

    % Floquet space dynamics 
    M = reshape(s(4*m+1:end), [m m]);                  %Monodromy matrix
    J = GNC.Control.PFSK.FloquetExponents;             %Floquet exponents of the reference trajectory
    F = GNC.Control.PFSK.FloquetDirections;            %Floquet exponents of the reference trajectory
    F = M*F*expm(-J*mod(t,T));                         %Instanteneous Floquet basis
    A = STM_dynamics(t, F, GNC);                       %Floquet and co-state dynamics

    % Final dynamics 
    ds = [ds; A*s(1:2*m)];
end

% Two-boundary value problem
function [dB] = bounds(x, x2, alpha0)
    %Constants 
    m = 6;          %Phase space dimension 

    %Initial state constraints
    dB(1:m) = x(1:m)-alpha0; 

    %Final state constraints
    dB(m+1:m+2) = [x2(m+1); -x2(m+2)]; 
    dB(m+3:m+4) = [1; -1];
    dB(11:60) = zeros(50,1);
end

% Variational dynamics
function [A] = STM_dynamics(t, F, GNC)
    % Constants 
    m = 6;          % Dimension of the STM 
    delta = 1e6;    % Saturation coefficient

    J = GNC.Control.PFSK.FloquetExponents;             %Floquet exponents of the reference trajectory
    lambda = GNC.Control.PFSK.InitialPrimer;           %Floquet modes of the reference trajectory
    cost_function = GNC.Control.PFSK.CostFunction;     %Cost function to minimize

    switch (cost_function)
        case 'L1'
            Tmax = GNC.Control.PFSK.MaxThrust;         %Maximum available thrust
        otherwise
            Tmax = 0;                                  %Maximum available thrust
    end

    % Compute the control law 
    B = [zeros(3,3); eye(3)];           %Control input matrix in the synodic frame

    %Compute the projected control matrix in the Floquet space 
    V = F^(-1)*B;
    p = -V.'*expm(-J.'*t)*lambda;       %Primer vector 

    %Switch depending on the cot function to minimize
    switch (cost_function)
        case 'L1' 
            %Final control vector
            Gamma = -V*V.'*Tmax/(1+exp(-delta*(norm(p)-1)));

         case 'L2'
            %Final control vector
            Gamma = -V*V.';          

         otherwise
            error('No valid vector norm to be minimized was selected');
    end

    % Compute the variational equations 
    A = [J Gamma; zeros(m,m) -J.']; 
end

% function [Sc, u, state] = PFSK_wrapper(mu, T, tf, s0, constraint, cost_function, Tmax)
%     %Constants 
%     m = 6;       %Phase space dimension
%     
%     %Sanity check on initial conditions dimension 
%     if (size(s0,1) == 1)
%         s0 = s0.';
%     end
%     
%     %Integration tolerance and setup 
%     options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
%     dt = 1e-3;                                                  %Integration time step  
%     
%     %Initial conditions and integration
%     Phi = eye(m);                                               %Initial STM 
%     Phi = reshape(Phi, [m^2 1]);                                %Reshape the initial STM
%     s0 = [s0; Phi];                                             %Complete phase space + linear variational initial conditions
%     
%     tspan = 0:dt:tf;                                            %Integration time span
%     [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
% 
%     %Energy constraint
%     constraint_flag = constraint.Flag;                          %Constraint flag
%     if (constraint_flag)
%         Cref = constraint.JacobiReference;                      %Reference Jacobi Constant
%     end
% 
%     %Definition of the GNC structure
%     GNC.Algorithms.Guidance = '';                       %Guidance algorithm
%     GNC.Algorithms.Navigation = '';                     %Navigation algorithm
%     GNC.Algorithms.Control = 'PFSK';                    %Control algorithm
%     
%     GNC.Guidance.Dimension = 6;                         %Dimension of the guidance law
%     GNC.Control.Dimension = 3;                          %Dimension of the control law
%     
%     GNC.System.mu = mu;                                 %Systems's reduced gravitational parameter
%     GNC.Control.PFSK.CostFunction = cost_function;      %Cost function to optimize
%     GNC.Control.PFSK.MaxThrust = Tmax;                  %Maximum available thurst
% 
%     %Floquet analysis                                             
%     [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), 0:dt:T, s0, options);
% 
%     Monodromy = reshape(Sn(end,2*m+1:end), [m,m]);              %Monodromy matrix of the relative orbit
%     [F,J] = eig(Monodromy);                                     %Eigenspectrum of the monodromy matrix
%     J = diag(log(diag(J))/T);                                   %Floquet exponents
%     for j = 1:size(F,2)
%         F(:,j) = F(:,j)/J(j,j);                                 %Floquet basis initial conditions
%     end
% 
%     GNC.Control.PFSK.FloquetExponents = J;                      %Floquet exponents of the reference trajectory
%     GNC.Control.PFSK.FloquetDirections = F;                     %Floquet directions of the reference trajectory
% 
%     %Differential corrector setup 
%     GoOn = true;                                                %Convergence boolean
%     maxIter = 100;                                              %Maximum initial iterations
%     iter = 0;                                                   %Initial iteration
%     tol = 1e-10;                                                %Differential corrector tolerance
% 
%     %Main computation
%     while (GoOn && (iter < maxIter))
%         %Evaluate the final error 
%         M = reshape(Saux(end,2*m+1:end), [m m]);                    %Instantenous Monodromy matrix
%         P = F*expm(-J*mod(tf,T));                                   %Full Floquet projection matrix
%         E = M*P;                                                    %Instantenous Floquet projection matrix
% 
%         Sf = Saux(end,1:m).'+Saux(end,m+1:2*m).';                   %Final absolute location
%         alpha(:,1) = zeros(6,1);                                    %Initial Floquet coordinates
%         alpha(:,2) = E^(-1)*Saux(end,m+1:2*m).';                    %Final Floquet coordinates
%         if (constraint_flag)
%             C = jacobi_constant(mu, Sf);                            %Final Jacobi Constant
%             error = [alpha(:,1); alpha(:,2); C-Cref];               %Error vector
%         else
%             error = [alpha(:,1); alpha(:,2)];                       %Error vector
%         end
% 
%         if (iter ~= 0)
%             %Sensitivity analysis 
%             A = STM_dynamics(tf, E, GNC);
%             STM = expm(A*tf);
%     
%             %Compute the differential corrector step 
%             ds = -STM\error; 
%     
%             %Re-evaluate the initial conditions 
%             lambda(:,2) = lambda(:,2) + ds(m+1:end);            %Initial primer vector conditions
%             GNC.Control.PFSK.InitialPrimer = lambda(:,2);       %Initial primer vector
%         else
%             %Initial co-state guess 
%             lambda(:,2) = expm(J.'*tf)*alpha(:,2); 
%             GNC.Control.PFSK.InitialPrimer = lambda(:,2);       %Initial primer vector
%         end
%         norm(error)
% 
%         %Re-integration of the trajectory
%         [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan, s0, options);
% 
%         %Convergence analysis
%         if (norm(error) < tol)
%             GoOn = false; 
%         else
%             iter = iter+1; 
%         end
%     end
% 
%     %Re-assembly of the control vector 
%     [~, ~, u] = GNCc_handler(GNC, Saux(:,1:m), Saux(:,m+1:end), tspan);
% 
%     %Outputs 
%     Sc = Saux(:,1:2*m);                     %Final trajectory output
%     state.State = ~GoOn;                    %Final convergence
%     state.Error = norm(error);              %Final error
%     state.Iter = iter;                      %Final iteration
% end

% function [Sc, u, state] = PFSK_wrapper(mu, T, tf, s0, constraint, cost_function, Tmax)
%     %Constants 
%     m = 6;       %Phase space dimension
%     
%     %Sanity check on initial conditions dimension 
%     if (size(s0,1) == 1)
%         s0 = s0.';
%     end
%     
%     %Integration tolerance and setup 
%     options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
%     dt = 1e-3;                                                  %Integration time step  
%     
%     %Initial conditions and integration
%     Phi = eye(m);                                               %Initial STM 
%     Phi = reshape(Phi, [m^2 1]);                                %Reshape the initial STM
%     s0 = [s0; Phi];                                             %Complete phase space + linear variational initial conditions
%     
%     tspan = 0:dt:tf;                                            %Integration time span
%     [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
% 
%     %Energy constraint
%     constraint_flag = constraint.Flag;                          %Constraint flag
%     if (constraint_flag)
%         Cref = constraint.JacobiReference;                      %Reference Jacobi Constant
%     end
% 
%     %Preallocation
%     u = zeros(3,length(tspan));                         %Control vector
%     S = zeros(size(Saux));                              %Close-loop trajectory
%     S(1,:) = s0; 
% 
%     %Definition of the GNC structure
%     GNC.Algorithms.Guidance = '';                       %Guidance algorithm
%     GNC.Algorithms.Navigation = '';                     %Navigation algorithm
%     GNC.Algorithms.Control = 'PFSK';                    %Control algorithm
%     
%     GNC.Guidance.Dimension = 6;                         %Dimension of the guidance law
%     GNC.Control.Dimension = 3;                          %Dimension of the control law
%     
%     GNC.System.mu = mu;                                 %Systems's reduced gravitational parameter
%     GNC.Control.PFSK.CostFunction = cost_function;      %Cost function to optimize
%     GNC.Control.PFSK.MaxThrust = Tmax;                  %Maximum available thurst
% 
%     %Main computation
%     for i = 1:length(tspan)-1
%         %Floquet analysis                                             
%         [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), 0:dt:T, [S(i,1:2*m).'; Phi], options);
% 
%         Monodromy = reshape(Sn(end,2*m+1:end), [m,m]);              %Monodromy matrix of the relative orbit
%         [F,J] = eig(Monodromy);                                     %Eigenspectrum of the monodromy matrix
%         J = diag(log(diag(J))/T);                                   %Floquet exponents
%         for j = 1:size(F,2)
%             F(:,j) = F(:,j)/J(j,j);                                 %Floquet basis initial conditions
%         end
%         P = F*expm(-J*mod(tspan(end)-tspan(i),T));                  %Full Floquet projection matrix
%     
%         GNC.Control.PFSK.FloquetExponents = J;                      %Floquet exponents of the reference trajectory
%         GNC.Control.PFSK.FloquetDirections = F;                     %Floquet directions of the reference trajectory
% 
%         %Evaluate the final error 
%         M = reshape(Saux(end,2*m+1:end), [m m]);                    %Instantenous Monodromy matrix
%         E = M*P;                                                    %Instantenous Floquet projection matrix
% 
%         Sf = Saux(end,1:m).'+Saux(end,m+1:2*m).';                   %Final absolute location
%         alpha = E^(-1)*Saux(end,m+1:2*m).';                         %Final Floquet coordinates
%         alpha(1)
%         if (constraint_flag)
%             C = jacobi_constant(mu, Sf);                            %Final Jacobi Constant
%             error = [alpha(1); C-Cref];                             %Error vector
%         else
%             error = alpha(1);                                       %Error vector
%         end
% 
%         %Stationkeeping constraints
%         if (norm(alpha(1:2)) ~= 0)
%             lambda(:,1) = [alpha(1:2)/norm(alpha(1:2)); zeros(4,1)];                      
%         else
%             lambda(:,1) = zeros(6,1);                       
%         end
% 
%         if (constraint_flag)
%             dC = jacobi_gradient(mu, Sf);                   %Gradient of the Jacobi Constant at the final time
%             lambda(:,1) = lambda(:,1) + E.'*dC;             %Final conditions on the primer vector
%         end
% 
%         %Evaluate the initial conditions 
%         lambda(:,2) = expm(J*(tf-tspan(i)))*lambda(:,1);    %Initial primer vector conditions
%         GNC.Control.PFSK.InitialPrimer = lambda(:,2);       %Initial primer vector
% 
%         %Re-integration of the trajectory
%         [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan(i:end), [S(i,1:2*m).'; Phi], options);
% 
%         %Re-assembly of the control vector 
%         [~, ~, du] = GNCc_handler(GNC, Saux(:,1:m), Saux(:,m+1:end), tspan);
%         u(:,i) = du(:,1);
% 
%         %Time-receding window 
%         S(i+1,:) = Saux(2,:);                 %Final trajectory
%     end
% 
%     u(:,i+1) = zeros(3,1);
% 
%     %Outputs 
%     Sc = S(:,1:2*m);                        %Final trajectory output
%     state.Error = norm(error);              %Final error
% end