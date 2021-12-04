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


function [Sc, u, state] = PFSK_wrapper(mu, T, tf, s0, constraint, cost_function, Tmax)
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
    P = F*expm(-J*mod(tf,T));                                   %Full Floquet projection matrix

    GNC.Control.PFSK.FloquetExponents = J;                      %Floquet exponents of the reference trajectory
    GNC.Control.PFSK.FloquetDirections = F;                     %Floquet directions of the reference trajectory

    %Differential corrector setup 
    GoOn = true;                                                %Convergence boolean
    maxIter = 100;                                              %Maximum initial iterations
    iter = 0;                                                   %Initial iteration

    %Main computation
    while (GoOn && (iter < maxIter))
        %Evaluate the final error 
        M = reshape(Saux(end,2*m+1:end), [m m]);                    %Instantenous Monodromy matrix
        E = M*P;                                                    %Instantenous Floquet projection matrix

        Sf = Saux(end,1:m).'+Saux(end,m+1:2*m).';                   %Final absolute location
        alpha = E^(-1)*Saux(end,m+1:2*m).';                         %Final Floquet coordinates
        alpha(1)
        if (constraint_flag)
            C = jacobi_constant(mu, Sf);                            %Final Jacobi Constant
            error = [alpha(1); C-Cref];                             %Error vector
        else
            error = alpha(1);                                       %Error vector
        end

        if (iter == 0)
            %Initial stationkeeping constraints guess 
            if (norm(error) ~= 0)
                lambda(:,1) = [alpha(1:2)/norm(alpha(1:2)); zeros(4,1)];                      
            else
                lambda(:,1) = zeros(6,1);                       
            end

            if (constraint_flag)
                dC = jacobi_gradient(mu, Sf);                   %Gradient of the Jacobi Constant at the final time
                lambda(:,1) = lambda(:,1) + E.'*dC;             %Final conditions on the primer vector
            end

            %Evaluate the initial conditions 
            lambda(:,2) = expm(J*tf)*lambda(:,1);               %Initial primer vector conditions
            GNC.Control.PFSK.InitialPrimer = lambda(:,2);       %Initial primer vector

        else
            %Sensitivity analysis 

            %Compute the differential corrector step 
            ds = -pinv(STM)*error; 
    
            %Re-evaluate the initial conditions 
            lambda(:,2) = lambda(:,2) + ds;                     %Initial primer vector conditions
            GNC.Control.PFSK.InitialPrimer = lambda(:,2);       %Initial primer vector
        end

        %Re-integration of the trajectory
        [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s, GNC), tspan, s0, options);

        %Convergence analysis
        if (norm(error) < tol)
            GoOn = false; 
        else
            iter = iter+1; 
        end
    end

    %Re-assembly of the control vector 
    [~, ~, u] = GNCc_handler(GNC, Saux(:,1:m), Saux(:,m+1:end), tspan);

    %Outputs 
    Sc = Saux(:,1:2*m);                     %Final trajectory output
    state.Error = norm(error);              %Final error
end

%% Auxiliary schemes
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