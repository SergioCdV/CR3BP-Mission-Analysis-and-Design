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
%         - vector Tmax, an initial guess for the maximum available thrust
%         - string solver, to select a T2BP solver of an algebraic one

% Output: - array Sc, the stationkeeping trajectory
%         - array u, containing the required control law
%         - structure state, with the corresponding differential corrector
%           figures of merit

% New versions: 

function [S, u, state] = PFSK_wrapper(mu, T, tf, s0, constraint, cost_function, Tmax, solver)
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
    S = zeros(size(Saux));                              %Close-loop trajectory
    S(1,:) = s0; 

    %Definition of the GNC structure
    GNC.Algorithms.Guidance = '';                       %Guidance algorithm
    GNC.Algorithms.Navigation = '';                     %Navigation algorithm
    GNC.Algorithms.Control = 'PFSK';                    %Control algorithm
    
    GNC.Guidance.Dimension = 6;                         %Dimension of the guidance law
    GNC.Control.Dimension = 3;                          %Dimension of the control law
    
    GNC.System.mu = mu;                                 %Systems's reduced gravitational parameter
    GNC.Control.PFSK.CostFunction = cost_function;      %Cost function to optimize
    GNC.Control.PFSK.MaxThrust = Tmax;                  %Maximum available thrust

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

    for i = 1:length(tspan)-1
        %Compute the initial guess        
        M = reshape(Saux(end,2*m+1:end), [m m]);                    %Instantenous Monodromy matrix
        P = F*expm(-J*mod(tspan(i),T));                             %Full Floquet projection matrix
        E = M*P;                                                    %Instantenous Floquet projection matrix
        alpha(:,1) = E^(-1)*Saux(1,m+1:2*m).';                      %Initial Floquet coordinates
        P = F*expm(-J*mod(tf-tspan(i),T));                          %Full Floquet projection matrix
        E = M*P;  
        alpha(:,2) = E^(-1)*Saux(end,m+1:2*m).';                    %Final Floquet coordinates
        lambda(:,2) = ones(6,1);                                    %Final co-state guess 
        lambda(:,1) = expm(J.'*(tf-tspan(i)))*lambda(:,2);          %Initial co-state guess
        GNC.Control.PFSK.InitialPrimer = lambda(:,2);               %Initial primer vector

        switch (solver)
            case 'BVP4C'
                %Initial guess for the numerical method
                nsteps = 10*length(0:dt:tf-tspan(i));
                solinit.x = linspace(0, tf-tspan(i), nsteps);
                solinit.y = repmat([alpha(:,1); lambda(:,1)], 1, nsteps);
            
                %Set optimizer options
                options = bvpset('RelTol', 1e-10, 'AbsTol', 1e-10*ones(12,1), 'Nmax', 2000);
            
                %Solve the problem 
                sol = bvp4c(@(t,s)dynamics(t, s, GNC), @(x,x2)bounds(x, x2, alpha(:,1)), solinit, options);
    
                %New initial primer vector
                GNC.Control.PFSK.InitialPrimer = sol.y(m+1:2*m,1);    

            case 'Newton'
                %Solve the dual minimum time/fuel problem usign a Newton method
                initial_guess = [alpha(:,1); lambda(:,1); tf; Tmax];
                sol = fsolve(@(t,s)pfskivp(t, s, alpha, GNC), initial_guess, options);

                %Solution setup
                GNC.Control.PFSK.InitialPrimer = sol(1:m);       %New initial primer vector   
                GNC.Control.PFSK.MaxThrust = sol(end);           %Maximum available thrust
                tspan = 0:dt:sol(end-1);                         %New integration time span

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
        dalpha(:,i) = alpha(:,1);
    end

    %Final step
    u(:,i+1) = zeros(3,1);

    %Outputs 
    S = S(:,1:2*m);            %Final trajectory output
    state.Error = dalpha;      %Final error
end

%% Auxiliary functions 
% Floquet stationkeeping IVP
function [e] = pfskivp(~, s, GNC, alpha0)
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
    [~, ~, u] = GNCc_handler(GNC, Saux(:,1:m), Saux(:,1:m), tspan);
    T = norm(u(:,end)); 

    %Compute the Hamiltonian of the problem evaluated at tf 
    alpha = Saux(end,1:m).';                        %Final Floquet coordinates
    lambda = Saux(end,m+1:2*m).';                   %Final co-state
    dS = dynamics(0, Saux(end,1:2*m).', GNC);       %Final dynamics vector field
    H = 1 + beta*T + dot(lambda, dS);               %Final problem Hamiltonian

    %Add the rest of the boundary conditions
    switch (problem)
        case 'Strict'
            e = [alpha(1); alpha(2:end)-alpha0(2:end); t-tf; H];                  
        case 'Minimization'
            e = [lambda(1)-1; lambda(2:end); t-tf; H];             
        case 'Rendezvous'
            e = [alpha; t-tf; H];                              
        otherwise
            error('No valid optimal control problem was selected')
    end
end

% Optimal problem dynamics 
function [ds] = dynamics(t, s, GNC)
    % Constants 
    m = 6;          %Phase space dimension 

    % Floquet space dynamics 
    F = GNC.Control.PFSK.FloquetDirections;            %Floquet exponents of the reference trajectory
    A = STM_dynamics(t, F, GNC);                       %Floquet and co-state dynamics

    % Final dynamics 
    ds = A*s(1:2*m);
end

% Two-boundary value problem constraints
function [dB] = bounds(x, x2, alpha0)
    %Constants 
    m = 6;          %Phase space dimension 

    %Initial state constraints
    dB(1:m) = x(1:m)-alpha0;

    %Final state constraints
    dB(m+1:2*m) = [x2(1); x2(2:6)-alpha0(2:end)];
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