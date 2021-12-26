%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 10/05/21
% File: FRTPS_control.m 
% Issue: 0 
% Validated: 26/12/21

%% Two impulsive, Floquet Rendezvous Target Mode Scheme %%
% This script contains the function to compute the control law by means of an FRTPS controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar TOF, the time of flight for the rendezvous condition
%         - vector s0, initial conditions of both the target and the
%           relative particle
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
%         - structure thruster_model, defining the thruster's error model
%         - scalar tol, the differential corrector scheme absolute

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [Sc, dV, state] = FRTPS_control(mu, TOF, s0, T, cost_function, sd, two_impulsive, penalties, target_points, thruster_model, tol)
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
    Rr = thruster_model.Rotation;        %Thrust error linear misalignement metric
    
    %Sanity check on the selected times indeces
    if (times(end) > size(Sn,1))
        warning('Final target point was located beyond the rendezvous time')
        times(end) = size(Sn,1); 
    end
    
    if (times(1) == 0)
        times(1) = 1; 
    end

    %Floquet analysis
    [~, Snaux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), 0:dt:T, s0, options);

    [E, lambda] = eig(reshape(Snaux(end,2*m+1:end), [m m]));
    J = diag(log(diag(lambda))/T);
    for i = 1:size(E,2)
        E(:,i) = E(:,i)/lambda(i,i);
    end
    F = inv(E);                         %Floquet change of basis

    u1 = [1; zeros(5,1)];               %Unstable direction vector
    B = u1.'*F*[zeros(m/2); Rr];        %Control input matrix

    %Initial error along the unstable manifold
    ns = u1.'*F*ns;
    
    %Implementation 
    while ((GoOn) && (iter < maxIter))
        %Compute the complete STM
        STM = reshape(S(end,13:end), [m m]);      %STM evaluated at time tf

        %Propagate the error 
        if (noise)
            gamma = 0;
            for i = 1:measurements 
                 gamma = gamma + M*exp(J(1,1)*tspan(times(i)));   %Noise state transition matrix
            end
            nSTM = sigma^2*B.'*gamma*B;                           %Accumulated state transition matrix
            nState = sigma*ns.'*gamma*B;                          %Accumulated noise vector
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
    dV = zeros(3,length(tspan));            %Impulses array
    dV(:,1) = (s0(10:12)-Sn(1,10:12).');    %Sum up each iteration contribution
    
    %Final impulse
    if (two_impulsive)
        dV(:,end) = -S(end,10:12).';        %Final impulse
        S(end,10:12) = zeros(1,3);          %Final conditions
    end

    Sc = S;                                 %Control trajectory 
    state.State = ~GoOn;                    %Convergence boolean
    state.Iterations = iter;                %Number of required iterations
    state.Error = norm(e);                  %Final error
end