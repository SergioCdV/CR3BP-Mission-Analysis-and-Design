%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/12/21
% File: CMC_guidance.m 
% Issue: 0 
% Validated: 06/12/21

%% Center Manifold Curve Guidance %%
% This script contains the function to compute the control law by means of the CMC guidance core.

% Inputs: - scalar mu, the reduced gravitational parameter of the system 
%         - scalar tf, the time of flight 
%         - scalar s0, the absolute target and chaser trajectory 
%         - scalar tol, the differential corrector tolerance

% Output: - array S, the converged guidance trajectory 
%         - structure state, information output about the differential
%           corrector process

% New versions: 

function [S, dV, state] = CMC_guidance(mu, T, tf, s0, tol)
    %Constants 
    m = 6;              %Relative phase space dimension
    dt = 1e-3;          %Time step
    tspan = 0:dt:tf;    %Integration time span

    %Integration tolerances 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

    %Integration of the relative chaser trajectory 
    s0c = [s0(1:m) s0(m+1:end)-s0(1:m)];                        %Initial relative chaser conditions
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0c, options);
    sf = Sn(end,m+1:end).';
            
    %Compute the initial state
    Phi = reshape(eye(m), [1 m^2]);             %Initial monodromy matrix
    epsilon = 1e-10;                            %Initial displacement along the center manifold
    s0 = [s0(1:m) zeros(1,m) Phi];              %Initial rendezvous conditions 
    [~, Sref] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), 0:dt:T, s0, options);

    Monodromy = reshape(Sref(end,2*m+1:end), [m m]);
    [E, Lambda] = eig(Monodromy);               %Eigenspectrum of the center manifold
    Ec = E(:,3);                                %Center manifold eigenvector
    Lambdac = Lambda(3,3);                      %Center manifold eigenvalue

    %Initial guess 
    rho = atan2(imag(Lambdac), real(Lambdac));  %Initial rotation number
    T0 = T;                                     %Initial stroboscopic time

    %Initial guess on the invariant curve 
    theta = pi/4; 
    s0(m+1:2*m) = epsilon*(real(Ec)*cos(theta).'-imag(Ec)*sin(theta).');           

    X = [s0(m+1:2*m).'; T0; rho];               %Complete initial guess

    Jref = jacobi_constant(mu, s0(1:m).');      %Reference Jacobi Constant

    [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
    [~, Scu] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), 0:dt:X(end-1), s0, options);

    %Differential corrector setup 
    iter = 1;                      %Initial iteration
    maxIter = 100;                 %Maximum number of iterations
    GoOn = true;                   %Convergence flag 

    while ((GoOn) && (iter < maxIter))
        %Invariance constraint
        R = exp(-1i*X(end));                            %Rotation operator
        utr = R*Scu(end,m+1:2*m).';                     %Rotated, stroboscopic state 
        dSc = utr-Scu(1,m+1:2*m).';                     %Phase invariance constraint 

        %Error analysis 
        Sf = Saux(1,1:m).'+Saux(1,m+1:2*m).';           %Initial current state
        J = jacobi_constant(mu, Sf);                    %Current Jacobi Constant
        dS = Saux(end,7:9).'-sf(1:3);                   %Final relative state difference
        error = [dSc; dS];                      %Error vector

        %Sensitivity analysis 
        Phi = reshape(Scu(end,2*m+1:end), [m m]);                                   %STM evaluated at the stroboscopic time
        M = kron(R,eye(m))*Phi-eye(m);                                              %STM matrix for the invariance constraint
        F = nlr_model(mu, true, false, false, 'Encke', 0, Scu(end,1:2*m).');        %Derivative with respect to the stroboscopic time
        F = R*F;                                                                    %Derivative with respect to the stroboscopic time
        dR = -1i*R*Scu(end,m+1:2*m).';                                              %Derivative with respect to the rotation number
        dJ = jacobi_gradient(mu, Sf);                                               %Jacobi Gradient
        STM = reshape(Saux(end,2*m+1:end), [m m]);                                  %State Transition Matrix
        A = [M F(7:12) dR; STM(1:3,:) zeros(3,2)];                                  %Sensitivity analysis

        %Newton-Rhapson update
        ds = -pinv(A)*error;            %Updat step
        X = X + ds;                     %New initial conditions
        s0(7:12) = X(1:m);              %New initial physical conditions
        norm(error)

        %Re-integration of the trajectory
        [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
        [~, Scu] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), 0:dt:X(end-1), s0, options);

        %Convergence analysis
        if (norm(error) < tol)
            GoOn = false; 
        else
            iter = iter + 1; 
        end
    end

    %Final output 
    dV(:,1) = Saux(1,10:12)-Sn(1,10:12);         %Initial impulse
    dV(:,2) = Saux(end,10:12)-Sn(end,10:12);     %Final impulse

    S = Saux;                   %Final trajectory
    state.State = ~GoOn;        %Final convergence flag 
    state.Iter = iter;          %Final iteration 
    state.Error = norm(error);  %Final error
end