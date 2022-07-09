%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date:  09/07/22
% File: FCG_control.m 
% Issue: 0 
% Validated: 09/07/22

%% Floquet Center Manifold Guidance Control %%
% This script contains the function to compute a guidance trajectory using
% the center manifold expansion of the relative dynamics

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar tf, the final time of flight
%         - vector s0, initial conditions of both the target and the chaser
%         - scalar tol, the differential corrector scheme tolerance for the
%           constrained maneuverç
%         - string restriction, specifying the type of required maneuver 

% Output: - array S, the rendezvous relative trajectory
%         - array dV, containing the required impulse located at the best
%           instant
%         - structure state, containing the results of the differential
%           corrector

% New versions: 

function [S, dV, state] = FCG_control(mu, tf, s0, tol, constraint, restriction)
    %Constants 
    m = 6;       % Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
    
    %Select the restriction level of the CAM 
    T = constraint.Period;                                      %Reference orbit period
    dt = 1e-3;                                                  %Integration time step  
    tspan = 0:dt:T;                                             %Integration time span
    
    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
        
    %Computation of the monodromy matrix
    s0 = [s0(1:m) s0(m+1:end)-s0(1:m)];                         %Initial relative chaser conditions
    Phi = eye(m);                                               %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                                %Reshape the initial STM
    s0 = [s0; Phi];                                             %Complete phase space + linear variational initial conditions

    [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 

    Monodromy = reshape(Saux(end,2*m+1:end), [m m]);            %Monodromy matrix
    [E, J] = eig(Monodromy);                                    %Eigenspectrum of the STM 

    for i = 1:size(J,2)
        E(:,i) = E(:,i)/J(i,i);
    end
    
    %Reference Jacobi constant
    Jref = jacobi_constant(mu, Saux(1,1:6).'+Saux(1,7:12).'); 
        
    %Differential corrector setup
    tspan = 0:dt:tf;                        %Integration time span
    maxIter = 100;                          %Maximum number of iterations

    %Initial integration 
    [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 
        
    %Differential corrector setup
    GoOn = true;                            %Convergence boolean
    iter = 1;                               %Initial iteration
        
    while ((GoOn) && (iter < maxIter))
        %Error vector
        error = [Saux(1,7:12).'; Saux(end,7:9).']; 

        Phi = reshape(Saux(end,2*m+1:end), [m m]);                      %State transition matrix at the final instant

        % Sensibility matrix for the center manifold restriction
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

        % Sensibility matrix for the rendezvous restriction
        C = Phi(1:3,:)*[E(:,3:end) [zeros(3,3); eye(3)]];
        STM = [STM; C];
        
        %Energy constraint 
        J = jacobi_constant(mu, s0(1:6).'+s0(7:12).');
        dJ = jacobi_gradient(mu, s0(1:6).'+s0(7:12).');
        JSTM = [zeros(1,size(STM,2)-3) dJ(4:6).'*STM(4:6,4:6)];

        %Complete sensibility analysis
        A = [STM; JSTM];                         %Sensibility matrix
        b = [error; J-Jref];                     %Final error analysis

        %Update the initial conditions
        ds = -pinv(A)*b;                         %Needed maneuver
        ds = ds(end-2:end);                      %Velocity change
        s0(10:12) = s0(10:12)+ds;                %Update initial conditions with the velocity impulse
     
        %Re-integrate the trajectory
        [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
        
        %Convergence analysis for the constrained case 
        if (norm(error) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end
     
    %Final output 
    dV(:,1) = Saux(1,10:12)-Sn(1,10:12);         %Initial impulse
    dV(:,2) = Saux(end,10:12)-Sn(end,10:12);     %Final impulse

    S = Saux;                                    %Final trajectory
    state.State = ~GoOn;                         %Final convergence flag 
    state.Iter = iter;                           %Final iteration 
    state.Error = norm(error);                   %Final error
end