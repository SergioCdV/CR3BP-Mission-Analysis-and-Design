%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/11/21
% File: MFSK_control.m 
% Issue: 0 
% Validated: 24/11/21

%% Modified Floquet Stationkeeping %%
% This script contains the function to compute the control law by means of an MFSK controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar T, the target orbit reference period
%         - vector s0, initial conditions of the target spacecraft
%         - scalar tol, the differential corrector scheme tolerance for the
%           constrained maneuver
%         - structure constraint, specifying any constraint on the maneuver
%         - vector Sg, reference initial orbit conditions
%         - scalar Jref, the Jacobi Constante reference value to be imposed

% Output: - array Sc, the stationkeeping trajectory
%         - array dV, containing the required impulse
%         - structure state, with the corresponding differential corrector
%           figures of merit

% New versions: 

function [Sc, dV, state] = MFSK_control(mu, T, s0, tol, constraint, Sg, Jref)
    %Constants 
    m = 6;       %Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
    
    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    dt = 1e-3;                                                  %Integration time step  
    tspan = 0:dt:T;                                             %Integration time span
    
    %Initial conditions and integration
    Phi = eye(m);                                               %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                                %Reshape the initial STM
    s0 = [s0; Phi];                                             %Complete phase space + linear variational initial conditions
    
    [~, Sn] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);
    S = Sn;
        
    %Differential corrector setup
    maxIter = 100;                          %Maximum number of iterations
    GoOn = true;                            %Convergence boolean
    iter = 1;                               %Initial iteration

    Constraint = constraint.Method;         %Method to constrain energy 
    constraint_flag = constraint.Flag;      %Constraint flag

    dV = zeros(11,1);                       %Initial iteration of the velocity impulse
        
    while ((GoOn) && (iter < maxIter))
        %Compute the Floquet modes at each time instant 
        Monodromy = reshape(S(end,m+1:end), [m m]);               %State transition matrix at each instant 
        [E, Lambda] = eig(Monodromy);                             %Eigenspectrum of the STM/Floquet basis 
        for j = 1:size(E,2)
            E(:,j) = E(:,j)/Lambda(j,j);
        end
            
        %Compute the maneuver
        error = s0(1:m)-Sg;                                       %State error vector
        STM = [E(:,2:end) -[zeros(3,3); eye(3)]];                 %Linear application
            
        %Energy constraint 
        if (constraint_flag)
            J = jacobi_constant(mu, s0(1:6));       %Actual Jacobi Constant
            switch (Constraint)
                case 'Impulse'
                    %Sensibility analysis
                    dJ = jacobi_gradient(mu, s0(1:6));
                    JSTM = [zeros(1,size(STM,2)-3) -dJ(4:6).'];
                    A = [STM; JSTM];
                    b = [error; J-Jref];

                    %Additional impulse
                    B = [zeros(3,3); eye(3)];
                    un = (E(:,2)-dot(E(:,2),E(:,1)));
                    dV2 = real(pinv(E^(-1)*B)*(1/2)*(J-Jref)/norm(s0(4:6))*un);  

                otherwise
                    %Sensibility analysis
                    dJ = jacobi_gradient(mu, s0(1:6));
                    JSTM = [zeros(1,size(STM,2)-3) -dJ(4:6).'];
                    A = [STM; JSTM];
                    b = [error; J-Jref];
            end
        else
            %Sensibility analysis
            A = STM;
            b = error;
        end
            
        %Compute the maneuver
        dV = pinv(A)*b;                   %Needed maneuver

        %Integrate the trajectory 
        switch (Constraint)
            case 'Impulse'
                s0(4:6) = s0(4:6)+dV(end-2:end)+dV2;               %Update initial conditions with the velocity impulse
            otherwise
                s0(4:6) = s0(4:6)+dV(end-2:end);                   %Update initial conditions with the velocity impulse
        end

        [~, S] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);
            
        %Convergence analysis for the constrained case 
        if (~constraint_flag)
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
    dV = (s0(4:6)-Sn(1,4:6).');
     
    %Integrate the SK trajectory
    s0 = Sn(1,1:m);                     
    [~, Sc] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, s0, options);

    %Differential corrector output
    state.State = ~GoOn;                    %Convergence boolean
    state.Iterations = iter;                %Number of required iterations
    state.Error = norm(error);              %Final error
end