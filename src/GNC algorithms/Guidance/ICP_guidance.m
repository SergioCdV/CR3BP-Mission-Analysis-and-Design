%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date:  09/07/22
% File: ICP_guidance.m 
% Issue: 0 
% Validated: 26/07/22

%% Impulsive Center Manifold Phasing %%
% This script contains the function to compute a phasing guidance trajectory using
% the center manifold expansion of the relative dynamics

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar tf, the final time of flight
%         - vector s0, initial conditions of both the target and the chaser
%         - scalar tol, the differential corrector scheme tolerance for the
%           constrained maneuver
%         - string restriction, specifying the type of required maneuver 

% Output: - array S, the rendezvous relative trajectory
%         - array dV, containing the required impulse located at the best
%           instant
%         - structure state, containing the results of the differential
%           corrector

% New versions: 

function [S, dV, state] = ICP_guidance(mu, L, gamma, T, dtheta, k, s0, tol, restriction)
    %Constants 
    m = 6;       % Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
    
    %Integration time span
    T = T*(1+dtheta/(2*pi*k));                                  %Needed relative period
    dt = 1e-3;                                                  %Time step  
    tspan = 0:dt:k*T;                                           %Integration time span

    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 

    %Reference Jacobi Constant
    Jref = jacobi_constant(mu, s0);

    %Computation of the monodromy matrix
    s0 = [s0(1:m); s0(m+1:end)-s0(1:m)];                        %Initial relative chaser conditions
    Phi = eye(m);                                               %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                                %Reshape the initial STM
    s0 = [s0; Phi];                                             %Complete phase space + linear variational initial conditions

    [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 
    sf = Saux(1,m+1:2*m);

    %Orbit parameters (frequencies)
    cn = legendre_coefficients(mu, L, gamma, 2);                %Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 %Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 %In-plane frequency
    wv  = sqrt(c2);                                             %Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 %Contraint on the planar amplitude

    %Compute the initial Lissajous 
    phi = atan2(sf(2)/kap, -sf(1));         %In-plane phase 
    Ax = -sf(1)/cos(phi);                   %In-plane amplitude 
    psi = -wp*k*T;                           %Out-of-plane phase 
    Az = sf(3)/sin(psi);                    %Out-of-plane amplitude

    %Seed trajectory
    s0(m+1:2*m) = zeros(1,m);               %Rendezvous conditions
    s0(7) = -Ax*cos(phi);                   %X relative coordinate
    s0(8) = kap*Ax*sin(phi);                %Y relative coordinate
    s0(9) = Az*sin(psi);                    %Z relative coordinate
    s0(10) = wp*Ax*sin(phi);                %Vx relative velocity
    s0(11) = kap*wp*Ax*cos(phi);            %Vy relative velocity
    s0(12) = wv*Az*cos(psi);                %Vz relative velocity

    %Initial integration
    [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options); 
    Sn = Saux; 
        
    %Differential corrector setup
    GoOn = true;                            %Convergence boolean
    maxIter = 50;                          %Maximum number of iterations
    iter = 1;                               %Initial iteration
        
    while ((GoOn) && (iter < maxIter))
        %K-th application of the stroboscopic map 
        Monodromy = reshape(Saux(2,2*m+1:end), [m m]);              %Monodromy matrix
        [E, J] = eig(Monodromy);                                    %Eigenspectrum of the initial monodromy matrix 
        for i = 1:size(J,2)
            E(:,i) = E(:,i)/J(i,i);
        end
    
        Phi = reshape(Saux(end,2*m+1:end), [m m]);                  %Monodromy matrix

        %Error vector
        error = [s0(7:12)]; 

        % Sensibility matrix for the center manifold restriction
        switch (restriction)
            case 'Mixed'
                STM = [E(:,2) E(:,3:end) -[zeros(3); eye(3)]];        %Linear application
            case 'Stable' 
                STM = [E(:,2) -[zeros(3); eye(3)]];                   %Linear application
            case 'Unstable'
                STM = [E(:,1) -[zeros(3); eye(3)]];                   %Linear application
            case 'Center'
                STM = [E(:,3:end) -[zeros(3); eye(3)]];               %Linear application
            otherwise
                error('No valid case was selected');
        end

        % Sensibility matrix for the rendezvous restriction
        C = [zeros(3,size(STM,2)-3) Phi(1:3,4:6)];
        
        %Energy constraint sensibility vector fields
        J = jacobi_constant(mu, Saux(1,1:6).'+Saux(1,7:12).');
        dJ = jacobi_gradient(mu, Saux(1,1:6).'+Saux(1,7:12).');
        JSTM = [zeros(1,size(STM,2)-3) dJ(4:6).'];

        %Complete sensibility analysis
        A = [STM; JSTM];                         %Sensibility matrix
        b = [error; Jref-J];                     %Final error analysis

        %Update the initial conditions
        ds = pinv(A)*b;                          %Needed maneuver
        dV = real(ds(end-2:end));                %Velocity change
        s0(10:12) = s0(10:12)+dV;                %Update initial conditions with the velocity impulse

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