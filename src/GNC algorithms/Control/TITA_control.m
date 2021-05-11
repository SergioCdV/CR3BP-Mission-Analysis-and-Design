%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 10/05/21
% File: TITA_control.m 
% Issue: 0 
% Validated: 10/05/21

%% Two impulsive, target approach control %%
% This script contains the function to compute the control law by means of an TITA controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar TOF, the time of flight for the rendezvous condition
%         - vector s0, initial conditions of both the target and the
%           relative particle
%         - scalar tol, the differential corrector scheme absolute
%           tolerance
%         - string cost_function, for both position, velocity and complete
%           rendezvous: 'Position', 'Velocity', 'State'
%         - vector sd, the desired state of the system at end of the
%           maneuver
%         - boolean two_impulsive, defining the two-impulsive/one-impulsive
%           strategy
%         - structure penalties, with the controller penalty matrices
%         - structure target_points, defining the target points for station
%           keeping
%         - structure thruster_model, defining the thruster's errro model

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [Sc, dV, state] = TITA_control(mu, TOF, s0, tol, cost_function, sd, two_impulsive, penalties, target_points, thruster_model)
    %Constants 
    m = 6;                               %Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
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
    
    %Controller characteristics 
    Qt = penalties.Q;                    %Penalty matrix on the state error
    R = penalties.R;                     %Penalty matrix on the control effort
    M = penalties.M;                     %Penalty matrix on stationkeeping
    
    Omegat = [zeros(3,3); eye(3)];       %Derivative of the state vector with respect to the impulse
        
    %Target points definition
    ns = 1e-6*ones(m,1);                 %Initial error state
    noise = target_points.Noise;         %Noise boolean to activate/desactivate error propagation
    times = target_points.Times;         %Times at which the maneuver shall be performed
    times = fix(times/dt);               %Position along the time span to impulse the spacecraft
    times = sort(times);                 %Sort the times at which the impulses are made
    
    measurements = length(times);        %Number of target points
    
    %Thruster model 
    sigma = thruster_model.Sigma;        %Thrust error linear scaling metric
    Rr = thruster_model.Rotation;        %Trhust error linear misalignement metric
    
    %Sanity check on the selected times indeces
    if (times(end) > size(Sn,1))
        warning('Final target point was located beyond the rendezvous time')
        times(end) = size(Sn,1); 
    end
    
    if (times(1) == 0)
        times(1) = 1; 
    end
    
    %Implementation 
    while ((GoOn) && (iter < maxIter))
        %Compute the complete STM
        STM = reshape(S(end,13:end), [m m]);      %STM evaluated at time tf

        %Propagate the error 
        if (noise)
            nSTM = zeros(m,m/2); 
            for i = 1:measurements 
                 dumbSTM = reshape(S(times(i),13:end), [m m]);              %Noise state transition matrix
                 nSTM = nSTM + dumbSTM.'*M*dumbSTM*[zeros(m/2); Rr];        %Accumulated state transition matrix
            end
            nSTM = sigma^2*[zeros(m,m/2) [zeros(m/2); Rr]]*nSTM;            %Accumulated state transition matrix
            nState = sigma*ns.'*nSTM;                                       %Accumulated noise vector
        end

        %Recompute initial conditions
        switch (cost_function)
            case 'Position' 
                xf = S(end,7:9);                %Final positon state
                Phi = STM(1:3,4:6);             %Position-velocity subSTM
                Q = Qt(1:3,1:3);                %Penalty matrix on the position error
                Omega = Omegat(4:6,:);          %Derivative of the state vector with respect to the impulse

                %Compute the STM 
                L = Omega.'*Phi.'*Q*Phi*Omega;  %Penalty on the state error 
                STM = R+L;                      %Sensibility matrix
                
            case 'Velocity'
                xf = S(end,10:12);              %Final positon state
                Phi = STM(4:6,4:6);             %Position-velocity subSTM
                Q = Qt(4:6,4:6);                %Penalty matrix on the position error
                Omega = Omegat(4:6,:);          %Derivative of the state vector with respect to the impulse

                %Compute the STM 
                L = Omega.'*Phi.'*Q*Phi*Omega;  %Penalty on the state error 
                STM = R+L;                      %Sensibility matrix

            case 'State' 
                xf = S(end,7:12);               %Final state
                Phi = STM(:,4:6);               %State-velocity transition matrix
                Q = Qt;                         %Penalty matrix on the complete state error
                Omega = Omegat(4:6,:);          %Derivative of the state vector with respect to the impulse

                %Compute the STM
                L = Omega.'*Phi.'*Q*Phi*Omega;  %Penalty on the state error 
                STM = R+L;                      %Sensibility matrix            

            otherwise
                error('No valid cost function was selected');
        end

        %Add some noise 
        if (noise)
            STM = STM + nSTM(4:6,1:3);                             %Noise state matrix
            e = (xf-sd)*Q*Phi*Omega + nState;                      %Error state (deviation from the rendezvous condition)
            dV(:,iter) = pinv(STM.')*e.';                          %Needed impulse
            s0(10:12) = s0(10:12)-dV(:,iter)+sigma*Rr*dV(:,iter);  %New initial conditions
            s0(7:12) = s0(7:12)+ns;                                %New noisy initial conditions
        else
            e = (xf-sd)*Q*Phi*Omega;                               %Error state (deviation from the rendezvous condition)
            dV(:,iter) = pinv(STM.')*e.';                          %Needed impulse
            s0(10:12) = s0(10:12)-dV(:,iter);                      %New initial conditions 
        end

        %Reintegrate the trajectory
        [~, S] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);

        %Convergence analysis 
        if (norm(e) < tol)
            GoOn = false;                        %Stop the method
        else
            iter = iter+1;                       %Update the iterations
        end
    end
    
    %Final initial impulse 
    dV1 = (s0(10:12)-Sn(1,10:12).');    %Sum up each iteration contribution
    
    %Final impulse
    if (two_impulsive)
        dV2 = -S(end,10:12).';          %Final impulse
        S(end,10:12) = zeros(1,3);      %Final conditions
    else
        dV2 = [];                       %No docking impulse
    end

    %Output
    dV = [dV1 dV2];                     %Impulses array
    Sc = S;                             %Control trajectory 
    state.State = ~GoOn;                %Convergence boolean
    state.Iterations = iter;            %Number of required iterations
    state.Error = norm(e);              %Final error

end