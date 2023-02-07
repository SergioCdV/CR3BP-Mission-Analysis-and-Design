%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 11/05/21
% File: FMSC_control.m 
% Issue: 0 
% Validated: 11/05/21

%% Floquet Mode Safe Control %%
% This script contains the function to compute the control law by means of an FMSC controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar TOC, the time of flight for the collision condition
%         - vector s0, initial conditions of both the target and the
%           relative particle
%         - structure setup, containing the definition of the algorithm

% Output: - vector tspan, the corresponding maneuver time span
%         - array Sc, the rendezvous relative trajectory
%         - array dV, containing the required impulse located at the best
%           instant

% New versions: 

function [tspan, S, dV] = FMSC_control(mu, TOC, s0, setup)
    % Compute the departure maneuver 
    [tspan, S, dV] = DEP_control(mu, TOC, s0, setup);

    % Compute the re-insertion maneuver and trajectory if allowed 
    if (setup.Reinsertion)
        [tspan, S, dV] = INS_control(mu, TOC, S, dV, setup);
    end
end

%% Auxiliary functions
% Planning and design of the re-insertion maneuver
function [tspan, S, dV] = INS_control(mu, TOC, S, dV, setup)
end 

% Planning and design of the departure maneuver
function [tspan, S, dV] = DEP_control(mu, TOC, s0, setup)
    % Constants 
    m = 6;                              % Phase space dimension
    STM_model = setup.STM;              % STM model to be used
    Ln = setup.LibrationID;             % Libration point identifier
    restriction = setup.Restriction;    % Dynamic structures to be used
    Tp = setup.ReferencePeriod;         % Reference orbit period
    
    % Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
        
    % Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      % Integration tolerances
    
    dt = 1e-3;                                                  % Integration time step  
    
    % Initial conditions and integration of the reference periodic relative orbit
    tspan = 0:dt:Tp;                                            % Integration time span

    Phi = eye(m);                                               % Initial STM 
    Phi = reshape(Phi, [m^2 1]);                                % Reshape the initial STM
    s0 = [s0; Phi];                                             % Complete phase space + linear variational initial conditions
    
    switch (STM_model)
        case 'Numerical'
            [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 

        case 'RLLM'
            % Relative and targe trajectories
            [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0(1:2*m), options); 

            % STM propagation 
            L = libration_points(mu);                       % System libration points
            gamma = L(end,Ln);                              % Characteristic distance of the libration point
            cn = legendre_coefficients(mu, Ln, gamma, 2);   % Legendre coefficient c_2 (equivalent to mu)
            c2 = cn(2);                                     % Legendre coefficient c_2 (equivalent to mu)
            Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];        % Potential linear term
            Omega = [0 2 0;-2 0 0; 0 0 0];                  % Coriolis force
            A = [zeros(m/2) eye(m/2); Sigma Omega];         % Jacobian matrix of the dynamics

            STM = zeros(length(tspan), m^2);
            for i = 1:length(tspan)
                STM = expm(A*tspan(i));
            end

            % Complete trajectory
            Sn = [Sn reshape(STM, [1 m^2])];

        otherwise
            error('No valid STM model was selected');
    end
    
    Jref = jacobi_constant(mu, (Sn(1,1:m) + Sn(1,m+1:2*m)).');     % Reference Jacobi constant
    Monodromy = reshape(Sn(end,2*m+1:end), [m m]);                 % State transition matrix at each instant 
    Phi = Monodromy;                                               % State transition matrix at the collision time
    [E, L] = eig(Monodromy);                                       % Eigenspectrum of the STM 
    E = E./diag(L);                                                % Floquet modes initial conditions

    % Compute the reference trajectory for the analysed TOF 
    tspan = 0:dt:TOC;                                           % Integration time span

    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 

    S = Sn;                                                     % Reference trajectory
        
    % Differential corrector setup
    maxIter = 100;                          % Maximum number of iterations
    GoOn = true;                            % Convergence boolean
    iter = 1;                               % Initial iteration
    tol = 1e-6;                             % Maneuver tolerance
        
    % Preallocation of the solution
    dV = zeros(3,length(tspan));            % Velocity impulse
                
    while ((GoOn) && (iter < maxIter))
        % Compute the maneuver
        e = s0(m+1:2*m);                                                % State error vector
        switch (restriction)
            case 'Mixed'
                STM = [E(:,1) E(:,3:end) -[zeros(3,3); eye(3)]];        % Linear application
            case 'Stable' 
                STM = [E(:,2) -[zeros(3,3); eye(3)]];                   % Linear application
            case 'Unstable'
                STM = [E(:,1) -[zeros(3,3); eye(3)]];                   % Linear application
            case 'Center'
                STM = [E(:,3:end) -[zeros(3,3); eye(3)]];               % Linear application
            otherwise
                error('No valid case was selected');
        end
        
        % Energy constraint 
        J = jacobi_constant(mu, s0(1:6)+s0(7:12));
        dJ = jacobi_gradient(mu, s0(1:6)+s0(7:12));
        JSTM = [zeros(1,size(STM,2)-3) dJ(4:6).'*Phi(4:6,4:6)];

        % Sensibility analysis
        A = [STM; JSTM];
        b = [e; Jref-J];
        
        % Compute the maneuver
        maneuver = real(pinv(A)*b);            % Needed maneuver
        dv = maneuver(end-2:end);              % Needed change in velocity

        % Integrate the trajectory 
        s0(10:12) = s0(10:12)+dv;              % Update initial conditions with the velocity impulse

        [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
        
        % Convergence analysis
        if (norm(dv) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end
        
    % Evaluate the final control law
    dV(:,1) = s0(10:12)-Sn(1,10:12).';
end

%% Previous versions
% First version
function [Sc, dV, tm] = FMSC_V1(mu, TOC, s0, tol, constraint, restriction)
    % Constants 
    m = 6;       % Phase space dimension
    
    % Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
    
    % Select the restriction level of the CAM 
    constrained = constraint.Constrained;                       %Select the type of maneuver
    lambda = constraint.SafeDistance;                           %Safety distance
    
    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    
    dt = 1e-3;                                                  %Integration time step  
    tspan = 0:dt:TOC;                                           %Integration time span
    
    %Initial conditions and integration
    Phi = eye(m);                                               %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                                %Reshape the initial STM
    s0 = [s0; Phi];                                             %Complete phase space + linear variational initial conditions
    
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 
    
    Jref = jacobi_constant(mu, Sn(1,1:6).'+Sn(1,7:12).');  
        
    %Differential corrector setup
    maxIter = 2;                               %Maximum number of iterations
    
    %Time horizon 
    time_horizon = length(tspan)-2;
    
    %Preallocation 
    dVf = zeros(3,time_horizon);                %Velocity impulse
        
    for i = 1:time_horizon
        %Integration set up 
        atime = tspan(i+1:end);                 %New time span
        s0 = Sn(i+1,:);                         %Initial conditions
        
        [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), atime, s0, options);
        
        %Differential corrector setup
        GoOn = true;                            %Convergence boolean
        iter = 1;                               %Initial iteration
        
        while ((GoOn) && (iter < maxIter))
            %Compute the Floquet modes at each time instant 
            if (constrained)
                Monodromy = reshape(Sn(end,13:end), [m m]);              %State transition matrix at each instant 
            else
                Monodromy = reshape(Sn(i,13:end), [m m]);                %State transition matrix at each instant
            end

            [E, L] = eig(Monodromy);                                     %Eigenspectrum of the STM 
            for j = 1:size(E,2)
                E(:,j) = E(:,j)/L(j,j);                                  %Compute the Floquet Modes
            end

            Phi = reshape(Sn(end,13:end), [m m]);                        %State transition matrix at the collision time
            
            %Compute the maneuver
            if (constrained)
                error = lambda*(E(1:3,1)+sum(E(1:3,3:end),2))-S(end,7:9).';     %Safety constraint
                STM = Monodromy(1:3,4:6);                                       %Correction matrix
            else
                error = s0(7:12).';                                             %State error vector
                switch (restriction)
                    case 'Mixed'
                        STM = [E(:,1) E(:,3:end) -[zeros(3,3); eye(3)]];        %Linear application
                    case 'Stable' 
                        STM = [E(:,2) -[zeros(3,3); eye(3)]];                   %Linear application
                    case 'Unstable'
                        STM = [E(:,1) -[zeros(3,3); eye(3)]];                   %Linear application
                    case 'Center'
                        STM = [E(:,3:end) -[zeros(3,3); eye(3)]];               %Linear application
                    otherwise
                        error('No valid case was selected');
                end
            end
            
            %Energy constraint 
            if (constraint.Energy)
                J = jacobi_constant(mu, s0(1:6).'+s0(7:12).');
                dJ = jacobi_gradient(mu, s0(1:6).'+s0(7:12).');
                JSTM = [zeros(1,size(STM,2)-3) dJ(4:6).'*Phi(4:6,4:6)];

                %Sensibility analysis
                A = [STM; JSTM];
                b = [error; Jref-J];
            else
                %Sensibility analysis
                A = STM;
                b = error;
            end
            
            %Compute the maneuver
            maneuver = pinv(A)*b;                  %Needed maneuver
            if (constrained)
                dV = real(maneuver);               %Needed change in velocity
            else
                dV = real(maneuver(end-2:end));    %Needed change in velocity
            end

            %Integrate the trajectory 
            s0(10:12) = s0(10:12)+dV.';            %Update initial conditions with the velocity impulse
            [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), atime, s0, options);
            
            %Convergence analysis for the constrained case 
            if (~constraint.Energy)
                GoOn = false;
            else
                if (norm(error) < tol)
                    GoOn = false;
                else
                    iter = iter+1;
                end
            end
        end
        
        %Evaluate the cost function
        dVf(:,i) = (s0(10:12)-Sn(i+1,10:12)).';
    end
     
    %Integrate the CAM trajectory
    [~, best] = sort(dot(dVf,dVf,1));                 %Select the minimum L2 norm solution over the time span
    index = 1; 
    GoOn = true;
    while (index <= size(dVf,2) && GoOn)
        if (norm(dVf(best(index))) ~= 0)
            GoOn = false;
        else
            index = index + 1;
        end
    end
    tm = tspan(best(index));                          %Time to perform the maneuver since detection
    atime = 0:dt:TOC-tm;                              %CAM integration time
    
    s0 = Sn(best(index),1:12);
    s0(10:12) = s0(10:12)+dVf(:,best(index)).';       %Update initial conditions with the velocity change
    
    dV = zeros(3,size(dVf,2));
    dV(:,best(index)) = dVf(:,best(index));
    [~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), atime, s0, options);
end