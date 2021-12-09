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
%         - string problem, the optimal problem to be solved
%         - scalar beta, the weight of the minimum fuel problem 
%         - string cost_function, to optimize on the l1 or l2 norm of the
%           control vector
%         - vector Tmax, an initial guess for the maximum available thrust
%         - string solver, to select a T2BP solver of an algebraic one

% Output: - array Sc, the stationkeeping trajectory
%         - array u, containing the required control law
%         - structure state, with the corresponding differential corrector
%           figures of merit

% New versions: 

function [S, u, state] = PFSK_wrapper(mu, T, tf, s0, constraint, problem, beta, cost_function, Tmax, solver)
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

    %Preallocation
    u = zeros(3,length(tspan));                         %Control vector
    dalpha = zeros(6,size(Saux,1));                     %Time history of the Floquet coordinates
    S = zeros(size(Saux));                              %Close-loop trajectory
    S(1,:) = s0;                                        %Initial conditions

    %Definition of the GNC structure
    GNC.Algorithms.Guidance = '';                       %Guidance algorithm
    GNC.Algorithms.Navigation = '';                     %Navigation algorithm
    GNC.Algorithms.Control = 'PFSK';                    %Control algorithm  
    GNC.Guidance.Dimension = 6;                         %Dimension of the guidance law
    GNC.Control.Dimension = 3;                          %Dimension of the control law
    GNC.System.mu = mu;                                 %Systems's reduced gravitational parameter
    GNC.Control.PFSK.CostFunction = cost_function;      %Cost function to optimize
    GNC.Control.PFSK.MaxThrust = Tmax;                  %Maximum available thrust
    GNC.Control.PFSK.Problem = problem;                 %Optimal problem to be solved
    GNC.Control.PFSK.Beta = beta;                       %Weigth of the minimum fuel problem
    GNC.Control.PFSK.Period = T;                        %Orbital period of the target orbit

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

    %Receding window
    for i = 1:length(tspan)-1
        %Compute the initial guess        
        M = reshape(Saux(end,2*m+1:end), [m m]);                    %Instantenous Monodromy matrix
        P = F*expm(-J*mod(tspan(i),T));                             %Full Floquet projection matrix
        E = M*P;                                                    %Instantenous Floquet projection matrix
        alpha(:,1) = E^(-1)*Saux(1,m+1:2*m).';                      %Initial Floquet coordinates
        lambda(:,1) = ones(m,1);                                   %Initial co-state guess
        GNC.Control.PFSK.InitialPrimer = lambda(:,1);               %Initial primer vector

        switch (solver)
            case 'BVP4C'
                %Initial guess for the numerical method
                nsteps = 10*length(0:dt:tf-tspan(i));
                solinit.x = linspace(0, tf-tspan(i), nsteps);
                solinit.y = repmat([alpha(:,1); lambda(:,1)], 1, nsteps);
            
                %Set optimizer options
                options = bvpset('RelTol', 1e-10, 'AbsTol', 1e-10*ones(12,1), 'Nmax', 2000);
            
                %Solve the problem 
                sol = bvp4c(@(t,s)dynamics(t, s, GNC), @(x,x2)bounds(x, x2, alpha(:,1), GNC), solinit, options);
    
                %New initial primer vector
                GNC.Control.PFSK.InitialPrimer = sol.y(m+1:2*m,1);    

            case 'Newton'
                %Solve the dual minimum time/fuel problem usign a Newton method
                initial_guess = [lambda(:,1); tf; Tmax];
                sol = fsolve(@(s)pfskivp(s, GNC, alpha(:,1)), initial_guess);

                %Solution setup
                GNC.Control.PFSK.InitialPrimer = sol(1:m);       %New initial primer vector   
                Tmax = sol(end);                                 %New maximum needed thrust
                GNC.Control.PFSK.MaxThrust = Tmax;               %Maximum available thrust
                tf = sol(end-1);                                 %New final time
                tspan = 0:dt:tf;                                 %New integration time span

            otherwise
                error('No valid solver was selected');
        end
    
        %Re-integration of the trajectory
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
        [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan(i:end), s0, options);
    
        %Re-assembly of the control vector 
        if (i ~= length(tspan)-1)
            [~, ~, du] = GNCc_handler(GNC, Saux(:,1:m), Saux(:,m+1:end), tspan(i:end));
        else
            du = zeros(3,1);
        end
    
        %Time-step solution 
        s0 = Saux(2,:);             %New initial conditions
        S(i+1,:) = Saux(2,:);       %New time step of the trajectory
        u(:,i) = du(:,1);           %Applied control vector
        dalpha(:,i) = alpha(:,1);   %Time history of the Floquet coordinates
    end

    %Final step
    u(:,i+1) = zeros(3,1);

    %Outputs 
    S = S(:,1:2*m);            %Final trajectory output
    state.Error = dalpha;      %Final error
end

%% Auxiliary functions 
% Floquet stationkeeping IVP
function [e] = pfskivp(s, GNC, alpha0)
    %Variables of interest 
    m = 6;                      %Phase space dimension
    tf = s(end-1);              %Minimum time of flight
    Tmax = s(end);              %Maximum thrust
    lambda = s(1:m);            %Initial co-state
    s0 = [alpha0; lambda];      %Initial conditions for the IVP

    GNC.Control.PFSK.MaxThrust = Tmax;                          %Maximum available thrust set up
    problem = GNC.Control.PFSK.Problem;                         %Optimal problem to be solved
    beta = GNC.Control.PFSK.Beta;                               %Weigth of the minimum fuel problem

    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    dt = 1e-3;                                                  %Integration time step 
    tspan = 0:dt:tf;                                            %Integration time span

    %Prepare and solve the IVP 
    [t, Saux] = ode113(@(t,s)dynamics(t, s, GNC), tspan, s0, options);

    %Evaluate the control law along the trajectory 
    J = GNC.Control.PFSK.FloquetExponents;                      %Floquet exponents of the reference trajectory
    F = GNC.Control.PFSK.FloquetDirections;                     %Floquet directions of the reference trajectory
    cost_function = GNC.Control.PFSK.CostFunction;              %Cost function to optimize
    T = GNC.Control.PFSK.Period;                                %Orbital period of the target orbit

    Sn = [zeros(1,m) reshape(eye(m), [1 m^2])]; 
    u = PFSK_control(t(end), T, Sn, F, J, Saux(end,m+1:2*m).', cost_function, Tmax);
    Th = norm(u(:,end)); 

    %Compute the Hamiltonian of the problem evaluated at tf 
    alpha = Saux(end,1:m).';                        %Final Floquet coordinates
    lambda = Saux(end,m+1:2*m).';                   %Final co-state
    dS = dynamics(0, Saux(end,1:2*m).', GNC);       %Final dynamics vector field
    H = beta*Th + dot(lambda, dS(1:m));             %Final problem Hamiltonian

    %Add the rest of the boundary conditions
    switch (problem)
        case 'Strict'
            e = [alpha(1); alpha(2:end)-alpha0(2:end); t(end)-tf; H];                  
        case 'Minimization'
            e = [lambda(1)-1; lambda(2:end); t(end)-tf; H];             
        case 'Rendezvous'
            e = [alpha; t(end)-tf; H];                              
        otherwise
            error('No valid optimal control problem was selected')
    end
end

% Optimal problem dynamics 
function [ds] = dynamics(t, s, GNC)
    % Floquet space dynamics 
    F = GNC.Control.PFSK.FloquetDirections;            %Floquet exponents of the reference trajectory
    ds = STM_dynamics(t, s, F, GNC);                   %Floquet and co-state dynamics
end

% Two-boundary value problem constraints
function [e] = bounds(x, x2, alpha0, GNC)
    %Constants 
    m = 6;                                  %Phase space dimension 
    problem = GNC.Control.PFSK.Problem;     %Optimal problem to be solved

    %Initial boundary conditions 
    e(1:m) = x(1:m)-alpha0; 

    %Final boundary conditions
    alpha = x2(1:m);                        %Final Floquet coordinates
    lambda = x2(m+1:2*m);                   %Final co-state 
    switch (problem)
        case 'Strict'
            e(m+1:2*m) = [alpha(1); alpha(2:end)-alpha0(2:end)];                  
        case 'Minimization'
            e(m+1:2*m) = [lambda(1)-1; lambda(2:end)];             
        case 'Rendezvous'
            e(m+1:2*m) = alpha;                              
        otherwise
            error('No valid optimal control problem was selected')
    end
end

% Variational dynamics
function [ds] = STM_dynamics(~, s, F, GNC)
    % Constants 
    m = 6;                  % Dimension of the STM 
    delta = 1e6;            % Saturation coefficient
    alpha = s(1:m);         % Floquet coordinates
    lambda = s(m+1:2*m);    % Co-state

    J = GNC.Control.PFSK.FloquetExponents;             % Floquet exponents of the reference trajectory
    cost_function = GNC.Control.PFSK.CostFunction;     % Cost function to minimize

    switch (cost_function)
        case 'L1'
            Tmax = GNC.Control.PFSK.MaxThrust;         %Maximum available thrust
        otherwise
            Tmax = 0;                                  %Maximum available thrust
    end

    % Compute the control law 
    B = [zeros(3,3); eye(3)];           %Control input matrix in the synodic frame

    %Compute the projected control matrix in the Floquet space 
    V = F^(-1)*B;          % Control input matrix
    p = -V.'*lambda;       % Primer vector 

    %Switch depending on the cot function to minimize
    switch (cost_function)
        case 'L1' 
            %Final control vector
            if (norm(p) ~= 0)
                Gamma = -V*(Tmax/(1+exp(-delta*(norm(p)-1))))*(p/norm(p));
            else
                Gamma = zeros(m,1);
            end

         case 'L2'
            %Final control vector
            Gamma = -V*p;          

         otherwise
            error('No valid vector norm to be minimized was selected');
    end

    % Compute the variational equations 
    ds = [J*alpha+Gamma; -J.'*lambda]; 
end