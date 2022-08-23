%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/12/21
% File: CMLG_guidance.m 
% Issue: 0 
% Validated: 06/12/21

%% Constrained Center Manifold Lissajous Guidance %%
% This script contains the function to compute the control law by means of the CMC guidance core.

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - scalar L, the Lagrange point around which the spacecraft orbits
%         - scalar gamma, the characteristic distance of the Lagrangian
%           point
%         - scalar tf, the time of flight 
%         - structure constraint, to impose requirements of the impulses
%           (to align with the center and stable manifolds)
%         - scalar s0, the absolute target and chaser trajectory 
%         - scalar tol, the differential corrector tolerance

% Output: - array S, the converged guidance trajectory 
%         - structure state, information output about the differential
%           corrector process

% New versions: 

function [S, dV, state, Sref] = CMLG_guidance(mu, L, gamma, tf, constraint, s0, tol)
    %Constants 
    m = 6;              %Relative phase space dimension
    dt = 1e-3;          %Time step
    tspan = 0:dt:tf;    %Integration time span

    %Integration tolerances 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

    %Orbit parameters (frequencies)
    cn = legendre_coefficients(mu, L(1), gamma, 2);             %Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 %Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 %In-plane frequency
    wv  = sqrt(c2);                                             %Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 %Contraint on the planar amplitude

    if (length(L) > 1)
        A = [zeros(3) eye(3); Sigma Theta]^(-1);                % Scaling matrix
    else
        A = eye(m);         % No scaling needed
    end

    %Integration of the relative chaser trajectory 
    s0c = [s0(1:m) s0(m+1:end)-A*s0(1:m)];                      %Initial relative chaser conditions
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0c, options);
    sf = Sn(1,m+1:end).';

    %Compute the initial center manifold state
    Phi = reshape(eye(m), [1 m^2]);                             %Initial monodromy matrix
    s0 = [s0(1:m) zeros(1,m) Phi];                              %Initial rendezvous conditions 

    if (constraint.Flag)
        T = constraint.Period;                                  %Synodic period of the relative orbit

        [~, Sref] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), 0:dt:T, s0, options);
    
        Monodromy = reshape(Sref(end,2*m+1:end), [m m]);            %Monodromy matrix
        [E, Lambda] = eig(Monodromy);                               %Eigenspectrum of the center manifold
        Ec = E(:,1)/Lambda(1,1);                                    %Unstable Floquet vector
    end
                      
    %Compute the initial Lissajous 
    phi = atan2(sf(2)/kap, -sf(1));         %In-plane phase 
    Ax = -sf(1)/cos(phi);                   %In-plane amplitude 
    psi = -wp*tf;                           %Out-of-plane phase 
    Az = sf(3)/sin(psi);                    %Out-of-plane amplitude

    %Seed trajectory
    s0(m+1:2*m) = zeros(1,m);               %Rendezvous conditions
    s0(7) = -Ax*cos(phi);                   %X relative coordinate
    s0(8) = kap*Ax*sin(phi);                %Y relative coordinate
    s0(9) = Az*sin(psi);                    %Z relative coordinate
    s0(10) = wp*Ax*sin(phi);                %Vx relative velocity
    s0(11) = kap*wp*Ax*cos(phi);            %Vy relative velocity
    s0(12) = wv*Az*cos(psi);                %Vz relative velocity

    tspan2 = 0:dt:2*tspan(end);
    r(1,:) = -Ax*cos(wp*tspan2+phi);         %X relative coordinate
    r(2,:) = kap*Ax*sin(wp*tspan2+phi);      %Y relative coordinate
    r(3,:) = Az*sin(wv*tspan2+psi);          %Z relative coordinate
    v(1,:) = wp*Ax*sin(wp*tspan2+phi);       %Vx relative velocity
    v(2,:) = kap*wp*Ax*cos(wp*tspan2+phi);   %Vy relative velocity
    v(3,:) = wv*Az*cos(wv*tspan2+psi);       %Vz relative velocity 

    [~, Saux] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan2, s0c(1:6), options);
    Sref = [Saux [r; v].'];

    [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);

    %Differential corrector setup 
    iter = 1;                      %Initial iteration
    maxIter = 100;                 %Maximum number of iterations
    GoOn = true;                   %Convergence flag 

    while ((GoOn) && (iter < maxIter))
        %Error analysis 
        dS = Saux(end,7:9).';                   %Final relative state difference

        if (constraint.Flag)
            V = [zeros(3,1); Saux(1,10:12).'-sf(4:6)];  %Initial velocity impulse
            dv = dot(V,Ec)-norm(Ec)*norm(V);            %Initial impulse constraint along the center manifold
            error = [dS; Saux(1,7:9).'; dv];            %Error vector
        else
            error = [dS; Saux(1,7:9).'-sf(1:3)];        %Error vector
        end

        %Sensitivity analysis 
        STM = reshape(Saux(end,2*m+1:end), [m m]);                      %State Transition Matrix
        A = [zeros(3) STM(1:3,4:6); eye(3) zeros(3)];                   %Derivative with respect to the state

        if (constraint.Flag)
            C = [A; zeros(1,3) Ec(4:6).'-norm(Ec)*V(1:3).'/norm(V)];    %Complete sensitivity matrix
        else
            C = A;                                                      %Complete sensitivity matrix
        end

        %Newton-Rhapson update
        ds = -pinv(C)*error;                   %Update step
        s0(7:12) = s0(7:12)+ds(1:m).';         %New initial conditions

        %Re-integration of the trajectory
        [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);

        %Convergence analysis
        if (norm(error(1:end-1)) < tol)
            GoOn = false; 
        else
            iter = iter + 1; 
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