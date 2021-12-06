%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/12/21
% File: CML_guidance.m 
% Issue: 0 
% Validated: 06/12/21

%% Center Manifold Lissajous Guidance %%
% This script contains the function to compute the control law by means of the CML guidance core.

% Inputs: - scalar mu, the reduced gravitational parameter of the system 
%         - scalar L, the Lagrange point ID around which the spacecraft
%           shall orbit
%         - scalar gamma, the characteristic distance of the Lagrange point
%         - scalar tf, the time of flight 
%         - scalar s0, the absolute target and chaser trajectory 
%         - scalar tol, the differential corrector tolerance

% Output: - array S, the converged guidance trajectory 
%         - structure state, information output about the differential
%           corrector process

% New versions: 

function [S, dV, state] = CML_guidance(mu, L, gamma, tf, s0, tol)
    %Constants 
    m = 6;              %Relative phase space dimension
    dt = 1e-3;          %Time step
    tspan = 0:dt:tf;    %Integration time span

    %Integration tolerances 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

    %Orbit parameters (frequencies)
    cn = legendre_coefficients(mu, L, gamma, 2);                %Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 %Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 %In-plane frequency
    wv  = sqrt(c2);                                             %Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 %Contraint on the planar amplitude

    %Integration of the relative chaser trajectory 
    s0c = [s0(1:m) s0(m+1:end)-s0(1:m)];                        %Initial relative chaser conditions
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, s0c, options);
    sf = Sn(end,m+1:end).';
            
    %Compute the initial Lissajous 
    Ax = norm(sf(1:2))/sqrt(1+kap^2);                           %In-plane amplitude 
    Az = sf(3)/sin(wv*tf);                                      %Out-of-plane amplitude
    psi = 0;                                                    %Out-of-plane phase 
    phi = atan2(sf(2),-sf(1)*kap)-wp*tf;                        %In-plane phase 

    %Seed trajectory
    s0(m+1:end) = zeros(1,m);           %Rendezvous conditions
    s0(7) = -Ax*cos(phi);               %X relative coordinate
    s0(8) = kap*Ax*sin(phi);            %Y relative coordinate
    s0(9) = Az*sin(psi);                %Z relative coordinate
    s0(10) = wp*Ax*sin(phi);            %Vx relative velocity
    s0(11) = kap*wp*Ax*cos(phi);        %Vy relative velocity
    s0(12) = wv*Az*cos(psi);            %Vz relative velocity

    Phi = reshape(eye(m), [1 m^2]);     %Initial STM
    s0 = [s0 Phi];                      %Initial conditions
    [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);

    %Differential corrector setup 
    iter = 1;                      %Initial iteration
    maxIter = 100;                 %Maximum number of iterations
    GoOn = true;                   %Convergence flag 

    while ((GoOn) && (iter < maxIter))
        %Error analysis 
        dS = Saux(end,7:9).'-sf(1:3);                   %Final relative state difference
        error = dS;                                     %Error vector

        %Sensitivity analysis 
        STM = reshape(Saux(end,2*m+1:end), [m m]);      %State Transition Matrix
        A = STM(1:3,4:6);                               %Sensitivity analysis

        %Newton-Rhapson update
        ds = -pinv(A)*error;                            %Updat step
        s0(10:12) = s0(10:12)+ds.';                     %New initial conditions

        %Re-integration of the trajectory
        [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);

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