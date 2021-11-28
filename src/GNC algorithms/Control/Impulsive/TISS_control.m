%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 10/05/21
% File: TISS_control.m 
% Issue: 0 
% Validated: 10/05/21

%% Two impulsive, single shooting control %%
% This script contains the function to compute the control law by means of an TISS controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar TOF, the time of flight for the rendezvous condition
%         - vector s0, initial conditions of both the target and the
%           relative particle
%         - scalar tol, the differential corrector scheme absolute
%           tolerance
%         - string cost_function, for both position, velocity and complete
%           rendezvous: 'Position', 'Velocity', 'State
%         - boolean two_impulsive, true if two impulses are allowed (one for
%           targetting, the second for docking)

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [Sc, dV, state] = TISS_control(mu, TOF, s0, tol, cost_function, two_impulsive)
    %Constants 
    m = 6;                               %Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
    
    %Sanity check on the cost function and the allowed number of impulses
    if (~two_impulsive)
        cost_function = 'State';
    end
    
    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    
    dt = 1e-3;                           %Integration time step  
    tspan = 0:dt:TOF;                    %Integration time span
    
    %Differential corrector set up
    maxIter = 100;                       %Maximum number of iterations
    GoOn = true;                         %Convergence boolean 
    iter = 1;                            %Initial iteration 
    
    %Preallocation 
    dV = zeros(3,maxIter);               %Targeting impulse

    %Initial conditions and integration
    Phi = eye(m);                        %Initial STM 
    Phi = reshape(Phi, [m^2 1]);         %Reshape the initial STM
    s0 = [s0; Phi];                      %Complete phase space + linear variational initial conditions
    
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
    S = Sn;
    
    %First impulse: targeting 
    while ((GoOn) && (iter < maxIter))
        %Rearrange the STM 
        STM = reshape(S(end,13:end), [m m]);                            %STM evaluated at time tf
        
        %Compute the error 
        switch (cost_function)
            case 'Position'
                error = S(end,7:9).';                                   %Position error
                STM = STM(1:3,4:6);                                     %Position-velocity subSTM
            case 'Velocity'
                error = S(end,10:12).';                                 %Velocity error
                STM = STM(4:6,4:6);                                     %Velocity-velocity subSTM
            case 'State'
                error = S(end,7:12).';                                  %Complete rendezvous error
                STM = STM(:,4:6);                                       %Complete state-velocity subSTM
            otherwise
                error('No valid cost function was chosen');
        end

        %Required impulse
        dV(1:3,iter) = pinv(STM)*error;                               

        %New initial conditions
        s0(10:12) = s0(10:12)-dV(:,iter);   
        
        %New trajectory
        [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 

        %Convergence analysis
        if (norm(error) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end
    
    %Final initial impulse 
    dV = zeros(3,length(tspan));            %Impulses array
    dV(:,1) = (s0(10:12)-Sn(1,10:12).');    %Sum up each iteration contribution
    
    %Final impulse
    if (two_impulsive)
        dV(:,end) = -S(end,10:12).';        %Final impulse
        S(end,10:12) = zeros(1,3);          %Final conditions
    end

    %Output
    Sc = S;                                 %Control trajectory 
    state.State = ~GoOn;                    %Convergence boolean
    state.Iterations = iter;                %Number of required iterations
    state.Error = norm(error);              %Final error
end
