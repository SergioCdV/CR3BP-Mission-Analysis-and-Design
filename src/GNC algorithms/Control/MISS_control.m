%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 10/05/21
% File: MISS_control.m 
% Issue: 0 
% Validated: 10/05/21

%% Multi impulsive, single shooting control %%
% This script contains the function to compute the control law by means of an MISS controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar TOF, the time of flight for the rendezvous condition
%         - vector s0, initial conditions of both the target and the
%           relative particle
%         - scalar tol, the differential corrector scheme absolute
%           tolerance
%         - string cost_function, for both position, velocity and complete
%           rendezvous: 'Position', 'Velocity', 'State
%         - structure impulses, defining the impulses along the trajectory.
%           The docking impulse is compulsory

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [Sc, dV, state] = MISS_control(mu, TOF, s0, tol, cost_function, impulses)
    %Constants 
    m = 6;       %Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
        
    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    
    dt = 1e-3;                                         %Integration time step  
    tspan = 0:dt:TOF;                                  %Integration time span
    
    %Differential corrector set up
    maxIter = 100;                                     %Maximum number of iterations
    GoOn = true;                                       %Convergence boolean 
    iter = 1;                                          %Initial iteration 
    
    %Initial conditions and integration
    Phi = eye(m);                                      %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                       %Reshape the initial STM
    s0 = [s0; Phi];                                    %Complete phase space + linear variational initial conditions
    
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
    S = Sn; 
    
    %Preallocation of the impulses 
    dV = zeros(m/2*impulses.Number, maxIter);          %Vector of impulses
    W = impulses.Weights;                              %Weightening matrix
    
    times = impulses.Times;                            %Times at which the maneuver shall be performed
    times = fix(times/dt);                             %Position along the time span to impulse the spacecraft
    times = sort(times);                               %Sort the times at which the impulses are made
    
    %Sanity check on the selected times indeces
    if (times(end) >= size(Sn,1))
        warning('Redundant docking impulse eliminated');
        times(end) = size(Sn,1)-1; 
    end
    
    if (times(1) == 0)
        times(1) = 1; 
    end
    
    impulses = impulses.Number;                        %Number of required impulses
    
    while ((GoOn) && (iter < maxIter))    
        %Compute the error
        switch (cost_function)
            case 'Position' 
                error = S(end,7:9).';                             %Final positon state     
                STM = zeros(m/2,m/2*impulses);                    %Preallocation of the STM

                %Compute the STM 
                for i = 1:length(times)
                    STMf = reshape(S(end,13:end), [m m]);         %STM evaluated at time tf
                    STM0 = reshape(S(times(i),13:end), [m m]);    %STM evaluated at time ti
                    STMdt = STMf*STM0^(-1);                       %STM evaluated between times tf-ti
                    STM(:,1+3*(i-1):3*i) = STMdt(1:3,4:6);        %Position-velocity subSTM
                end

            case 'Velocity' 
                error = S(end,10:12).';                           %Final velocity state
                STM = zeros(m/2,m/2*impulses);                    %Preallocation of the STM

                %Compute the STM 
                for i = 1:length(times)
                    STMf = reshape(S(end,13:end), [m m]);         %STM evaluated at time tf
                    STM0 = reshape(S(times(i),13:end), [m m]);    %STM evaluated at time ti
                    STMdt = STMf*STM0^(-1);                       %STM evaluated between times tf-ti
                    STM(:,1+3*(i-1):3*i) = STMdt(4:6,4:6);        %Velocity-velocity subSTM
                end

            case 'State' 
                error = S(end,7:12).';                            %Final state error
                STM = zeros(m,m/2*impulses);                      %Preallocation of the STM

                %Compute the STM 
                for i = 1:length(times)
                    STMf = reshape(S(end,13:end), [m m]);         %STM evaluated at time tf
                    STM0 = reshape(S(times(i),13:end), [m m]);    %STM evaluated at time ti
                    STMdt = STMf*STM0^(-1);                       %STM evaluated between times tf-ti
                    STM(:,1+3*(i-1):3*i) = STMdt(:,4:6);          %State-velocity subSTM
                end

            otherwise
                error('No valid cost function was selected');
        end

        %Control law (impulses law)
        dV(:,iter) = -inv(W)*STM.'*(STM*W*STM.')^(-1)*error;                   

        %Reintegrate the trajectory
        for i = 1:length(times)               
            %Make the impulse
            S(times(i),10:12) = S(times(i),10:12) + dV(1+3*(i-1):3*i,iter).';

            %Integration time span
            if (i ~= length(times))
                Dt = tspan(times(i)):dt:tspan(times(i+1));
            else
                Dt = tspan(times(i)):dt:tspan(end);
            end

            %New trajectory
            [~, s] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), Dt, S(times(i),:), options); 

            %Update new initial conditions
            S(times(i):times(i)+size(s,1)-1,:) = s;
        end

        %Convergence analysis 
        if (norm(error) < tol)
            GoOn = false;                       %Stop the method
        else
            iter = iter+1;                      %Update the iterations
        end
    end
            
    %Output       
    dVf = dV(:,iter);                                   %Final impulses
    dV = zeros(3,length(tspan));                        %Impulses array
    dV(:,times) = reshape(dVf, [3 impulses]);           %Reshape the impulses array
    dV(:,1) = (S(1,10:12)-Sn(1,10:12)).';               %Initial impulse
    dV(:,end) = -S(end,10:12).';                        %Final impulse, always needed for docking
    Sc = S;                                             %Control trajectory 
    Sc(end,10:12) = zeros(1,3);                         %Nullify the final relative velocity
    state.State = ~GoOn;                                %Convergence boolean
    state.Iterations = iter;                            %Number of required iterations
    state.Error = norm(error);                          %Final error
end