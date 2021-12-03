%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/11/21
% File: PFSK_control.m 
% Issue: 0 
% Validated: 30/11/21

%% Primer Floquet Stationkeeping %%
% This script contains the function to compute the control law by means of an PFSK controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar T, the target orbit reference period
%         - vector s0, initial conditions of the target spacecraft
%         - scalar tol, the differential corrector scheme tolerance for the
%           constrained maneuver
%         - structure constraint, specifying any constraint on the maneuver
%         - string cost_function, to optimize on the l1 or l2 norm of the
%           control vector
%         - vector Tmax, the maximum available thrust

% Output: - array Sc, the stationkeeping trajectory
%         - array u, containing the required control law
%         - structure state, with the corresponding differential corrector
%           figures of merit

% New versions: 

function [Sc, u, state] = PFSK_control(mu, T, tf, s0, tol, constraint, cost_function, Tmax)
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
    
    tspan = 0:dt:T;                                             %Integration time span
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);

    %Floquet analysis                                           %Orbit-period time span               
    Monodromy = reshape(Sn(end,2*m+1:end), [m,m]);              %Monodromy matrix of the relative orbit
    [F,J] = eig(Monodromy);                                     %Eigenspectrum of the monodromy matrix
    J = diag((1/T)*log(diag(J)));                               %Floquet exponents
    for i = 1:size(F,2)
        F(:,i) = F(:,i)/J(i,i);                                 %Floquet basis initial conditions
    end

    tspan = 0:dt:tf;                                            %Integration time span
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
    S = Sn;
        
    %Differential corrector setup
    maxIter = 100;                                              %Maximum number of iterations
    GoOn = true;                                                %Convergence boolean
    iter = 1;                                                   %Initial iteration

    %Energy constraint
    constraint_flag = constraint.Flag;                          %Constraint flag
    if (constraint_flag)
        Jref = constraint.JacobiReference;                      %Reference Jacobi Constant
    end

    %Preallocation of the control vector 
    u = zeros(3,length(tspan));

    %Definition of the GNC structure
    GNC.Algorithms.Guidance = '';                       %Guidance algorithm
    GNC.Algorithms.Navigation = '';                     %Navigation algorithm
    GNC.Algorithms.Control = 'OPFSK';                   %Control algorithm
    
    GNC.Guidance.Dimension = 9;                         %Dimension of the guidance law
    GNC.Control.Dimension = 3;                          %Dimension of the control law
    
    GNC.System.mu = mu;                                 %Systems's reduced gravitational parameter
    GNC.Control.OPFSK.CostFunction = cost_function;     %Cost function to optimize
    GNC.Control.OPFSK.FloquetExponents = J;             %Floquet exponents of the reference trajectory
    GNC.Control.OPFSK.FloquetDirections = F;            %Floquet directions of the reference trajectory
    GNC.Control.OPFSK.MaxThrust = Tmax;                 %Maximum available thurst

    %Main computation
    while ((GoOn) && (iter < maxIter))
        %Evaluate the final error 
        Sf = S(end,1:m).'+S(end,m+1:2*m).';                 %Final absolute location
        M = reshape(S(end,2*m+1:end), [m m]);               %Instantenous Monodromy matrix
        E = M*F*expm(-J*mod(tf,T));                         %Instantenous Floquet projection matrix
        dC = jacobi_gradient(mu, Sf);                       %Gradient of the Jacobi Constant at the final time

        lambda(:,1) = [ones(1,1); zeros(5,1)];              %Stationkeeping constraints
        if (constraint_flag)
            lambda(:,1) = lambda(:,1) + E.'*dC;             %Final conditions on the primer vector
        end

        %Error analysis
        alpha = E^(-1)*S(end,m+1:2*m).';                    %Final Floquet coordinates
        error = alpha(1:2);                                 %Error vector
        if (constraint_flag)
            J = jacobi_constant(mu, Sf);                    %Final Jacobi Constant
            error = [error; J-Jref];
        end

        %Evaluate the initial conditions 
        lambda(:,2) = expm(J*tf)*lambda(:,1);               %Initial primer vector conditions
        GNC.Control.OPFSK.InitialPrimer = lambda(:,2);      %Initial primer vector

        %Re-integration of the trajectory
        [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan, s0, options);

        %Re-assembly of the control vector 
        [~, ~, du] = GNCc_handler(GNC, S(:,1:m), S(:,m+1:end), tspan);
        u(:,1:size(du,2)) = u(:,1:size(du,2)) + du;
            
        %Convergence analysis 
        if (norm(error) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end

    %Outputs 
    Sc = S(:,1:2*m);                        %Final trajectory output
    state.State = ~GoOn;                    %Convergence boolean
    state.Iterations = iter;                %Number of required iterations
    state.Error = norm(error);              %Final error
end