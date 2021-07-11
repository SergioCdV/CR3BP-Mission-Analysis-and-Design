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
%         - vector so, the relative state of the colliding object
%         - vector s0, initial conditions of both the target and the
%           relative particle
%         - matrix Q, a positive definite matrix projecting the quadratic
%           form of the error to the colliding object
%         - scalar tol, the differential corrector scheme tolerance for the
%           constrained maneuver
%         - structure constraint, specifying any constraint on the maneuver
%         - string resctriction, specifying the type of required maneuver 

% Output: - array Sc, the rendezvous relative trajectory
%         - array dVf, containing the required impulse located at the best
%           instant

% New versions: 

function [Sc, dVf, tm] = FMSC_control(mu, TOC, so, s0, Q, tol, constraint, restriction)
    %Constants 
    m = 6;       %Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
    
    %Select the restriction level of the CAM 
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
    
    %Preallocation 
    J = zeros(2,length(tspan)-1);                               %Cost function to analyze
    dVf = zeros(3,length(tspan)-1);                             %Velocity impulse
    
    %Differential corrector setup
    maxIter = 20;                                               %Maximum number of iterations
    
    %Time horizon 
    if (constrained)
        time_horizon = length(tspan);
    else
        time_horizon = 1;
    end
        
    for i = 1:time_horizon
        %Integration set up 
        atime = tspan(i:end);                   %New time span
        s0 = Sn(i,:);                           %Initial conditions
        
        [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), atime, s0, options);
        
        %Differential corrector setup
        GoOn = true;                            %Convergence boolean
        iter = 1;                               %Initial iteration
        
        while ((GoOn) && (iter < maxIter))
            %Compute the Floquet modes at each time instant 
            Monodromy = reshape(S(end,13:end), [m m]);                     %State transition matrix at each instant 
            Monodromy = Monodromy*reshape(S(1,13:end), [m m])^(-1);        %Relative STM

            [E, ~] = eig(Monodromy);                                       %Eigenspectrum of the STM 
            for j = 1:size(E,2)
                E(:,j) = E(:,j)/norm(E(:,j));                              %Normalize unit basis
            end
            
            STM = reshape(S(1,13:end), [m m]);
            [E2,~] = eig(STM);

            %Compute the maneuver
            if (constrained)
                error = lambda*(E(1:3,1)+sum(E(1:3,3:end),2))-S(end,7:9).';     %Safety constraint
                STM = Monodromy(1:3,4:6);                                       %Correction matrix
            else
                error = s0(7:12).';                                             %State error vector
                switch (restriction)
                    case 'Best'
                        STM = [E(:,1) E(:,3:end) -[zeros(3,3); eye(3)]];        %Linear application
                    case 'Worst'
                        STM = [E(:,1) -[zeros(3,3); eye(3)]];                   %Linear application
                    otherwise
                        error('No valid case was selected');
                end
            end

            %Compute the maneuver
            maneuver = pinv(STM)*error;            %Needed maneuver
            if (constrained)
                dV = real(maneuver);               %Needed change in velocity
            else
                dV = real(maneuver(end-2:end));    %Needed change in velocity
            end

            %Integrate the trajectory 
            s0(10:12) = s0(10:12)+real(dV).';  %Update initial conditions with the velocity impulse
            [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), atime, s0, options);
            
            %Convergence analysis for the constrained case 
            if (constrained)
                if (norm(error) < tol)
                    GoOn = false;
                else
                    iter = iter+1;
                end
            else
                GoOn = false;
            end
        end
        
        %Evaluate the cost function
        dVf(:,i) = (s0(10:12)-Sn(i,10:12)).';
        J(1,i) = (1/2)*(S(end,7:9)-so(1,1:3))*Q*(S(end,7:9)-so(1,1:3)).';
        J(2,i) = norm(dVf(:,i))+(1/2)*S(end,7:12)*S(end,7:12).';  
    end
     
    %Select the minimum dV maneuver maximizing the relative distance to the colliding object
    tol = 1e-10;
    best = 1; 
    for i = 1:size(J,2)
        if (J(1,i)-J(1,best) > tol)
            best = i;
        elseif (J(1,i)-J(1,best) < tol)
            %Do nothing, just the other extreme case
        else
            %Minimize the second cost function 
            if (J(2,i) < J(2,best))
                best = i;
            end
        end
    end

    %Integrate the CAM trajectory
    tm = tspan(best(end));                            %Time to perform the maneuver since detection
    atime = 0:dt:TOC-tm;                              %CAM integration time
    
    if (isempty(atime))
        atime = tspan;
        s0 = Sn(1,1:12); 
        s0(10:12) = s0+dVf(:,1).';                    %Update initial conditions with the velocity change
    else
        s0 = Sn(best(end),1:12);
        s0(10:12) = s0(10:12)+dVf(:,best(end)).';     %Update initial conditions with the velocity change
    end
    
    dVf = dVf(:,best(end));
    [~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), atime, s0, options);
end