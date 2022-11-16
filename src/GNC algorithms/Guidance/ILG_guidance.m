%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/12/21
% File: ILG_guidance.m 
% Issue: 0 
% Validated: 06/12/21

%% Constrained Center Manifold Impulsive Lissajous Guidance %%
% This script contains the function to compute the control law by means of the ILG guidance core

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

function [S, dV, state, Sref, SM] = ILG_guidance(mu, L, gamma, tf, constraint, s0, tol)
    % Constants 
    m = 6;                                                      % Relative phase space dimension
    dt = 1e-3;                                                  % Time step
    tspan = 0:dt:tf;                                            % Integration time span

    % Integration tolerances 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

    % Orbit parameters
    cn = legendre_coefficients(mu, L(1), gamma, 2);             % Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 % Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 % In-plane frequency
    wv  = sqrt(c2);                                             % Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 % Constraint on the planar amplitude

    % Transformation to the affine libration point
    if (length(L) > 1)
        Theta = [0 2 0; -2 0 0; 0 0 0];
        Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];
        SM = [zeros(3) eye(3); Sigma Theta]^(-1);               % Scaling matrix

        cn = legendre_coefficients(mu, L(2), gamma, 2);         % Legendre coefficient c_2 (equivalent to mu)
        c2 = cn(2);                                             % Legendre coefficient c_2 (equivalent to mu)
        Theta = [0 2 0; -2 0 0; 0 0 0];
        Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];
        SM = SM*[zeros(3) eye(3); Sigma Theta];                 % Scaling matrix
    else
        SM = eye(m);                                            % No scaling needed
    end
    
    % Initial conditions 
    sd = s0(m+1:2*m).'-SM*s0(1:m).';                            % Initial affine Lissajous phase-space vector
    s0 = [s0(1:m) s0(m+1:2*m)-s0(1:m)];                         % Initial relative chaser conditions
    Phi = reshape(eye(m), [1 m^2]);                             % Initial monodromy matrix
    s0 = [s0 Phi];                                              % Complete initial conditions

    % Initial rendezvous guess
    [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
    Sn = Saux;

    % Center manifold constraints
    if (constraint.Flag)
        T = constraint.Period;                                  % Synodic period of the relative orbit
        [~, Sref] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), 0:dt:T, s0, options);
    
        Monodromy = reshape(Sref(end,2*m+1:end), [m m]);        % Monodromy matrix
        [E, Lambda] = eig(Monodromy);                           % Eigenspectrum of the center manifold
        Ec = E(:,1)/Lambda(1,1);                                % Unstable Floquet vector
    end
                      
    % Compute the initial Lissajous trajectory seed
    phi = atan2(sd(2)/kap, -sd(1));                                 % In-plane phase 
    Ax = -sd(1)/cos(phi);                                           % In-plane amplitude 
    psi = -wp*tf;                                                   % Out-of-plane phase 
    Az = sd(3)/sin(psi);                                            % Out-of-plane amplitude

    % Lissajous curve guess
    tspan2 = 0:dt:2*tspan(end);
    [~, Slis] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan2, s0(1:6), options);

    r(1,:) = -Ax*cos(wp*tspan2+phi);         % X relative coordinate
    r(2,:) = kap*Ax*sin(wp*tspan2+phi);      % Y relative coordinate
    r(3,:) = Az*sin(wv*tspan2+psi);          % Z relative coordinate
    v(1,:) = wp*Ax*sin(wp*tspan2+phi);       % Vx relative velocity
    v(2,:) = kap*wp*Ax*cos(wp*tspan2+phi);   % Vy relative velocity
    v(3,:) = wv*Az*cos(wv*tspan2+psi);       % Vz relative velocity 

    Sref = [Slis [r; v].'];                       % Complete relative state phase-space trajectory

    % Lissajous guess
    Saux(:,7) = -Ax*cos(wp*tspan+phi).';          % X relative coordinate
    Saux(:,8) = kap*Ax*sin(wp*tspan+phi).';       % Y relative coordinate
    Saux(:,9) = Az*sin(wv*tspan+psi).';           % Z relative coordinate
    Saux(:,10) = wp*Ax*sin(wp*tspan+phi).';       % Vx relative velocity
    Saux(:,11) = kap*wp*Ax*cos(wp*tspan+phi).';   % Vy relative velocity
    Saux(:,12) = wv*Az*cos(wv*tspan+psi).';       % Vz relative velocity

    s0(m+1:2*m) = Saux(1,m+1:2*m);

    % Differential corrector setup 
    iter = 1;                                % Initial iteration
    maxIter = 100;                           % Maximum number of iterations
    GoOn = true;                             % Convergence flag 

    while ((GoOn) && (iter < maxIter))
        % Error analysis 
        if (constraint.Flag)
            V = [zeros(3,1); Saux(1,10:12).'-sd(4:6)];  % Initial velocity impulse
            dv = dot(V,Ec)-norm(Ec)*norm(V);            % Initial impulse constraint along the center manifold
        else
            dv = [];                                    % Null initial impulse constraint along the center manifold      
        end

        error = [Saux(end,7:9).'; Saux(1,7:9).'-sd(1:3); dv]; 

        % Sensitivity analysis 
        STM = reshape(Saux(end,2*m+1:end), [m m]);                      % State Transition Matrix
        A = [zeros(3) STM(1:3,4:6); eye(3) zeros(3)];                   % Sensibillity matrix

        if (constraint.Flag)
            C = [A; zeros(1,3) Ec(4:6).'-norm(Ec)*V(1:3).'/norm(V)];    % Complete sensibillity matrix
        else
            C = A;                                                      % Complete sensibillity matrix
        end
 
        % Newton-Rhapson update
        ds = -pinv(C)*error;                  % Update step
        s0(7:12) = s0(7:12)+ds(1:m).';        % New initial conditions

        % Re-integration of the trajectory
        [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);

        % Convergence analysis
        if (norm(error) < tol)
            GoOn = false; 
        else
            iter = iter + 1; 
        end
    end

    % Final output 
    dV = zeros(3,length(tspan));                 % Impulses array
    dV(:,1) = Saux(1,10:12)-Sn(1,10:12);         % Initial impulse
    dV(:,end) = -Saux(end,10:12);                % Final impulse

    S = Saux;                                    % Final trajectory
    S(end,10:12) = zeros(1,3);                   % Final null relative velocity
    state.State = ~GoOn;                         % Final convergence flag 
    state.Iter = iter;                           % Final iteration 
    state.Error = norm(error);                   % Final error
end