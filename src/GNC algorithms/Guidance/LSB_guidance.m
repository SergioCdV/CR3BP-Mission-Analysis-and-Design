%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date:  09/07/22
% File: ICP_guidance.m 
% Issue: 0 
% Validated: 26/07/22

%% Lissajous Shape-based Guidance %%
% This script contains the function to compute a phasing guidance trajectory using
% the center manifold expansion of the relative dynamics

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar L, the ID of the libration point around which the target
%           orbits
%         - scalar gamma, the characteristic distance of the libration
%           point
%         - vector s0, initial conditions of both the target and the chaser
%         - string method, the solver of the optimal control TVBP
%         - scalar tf, the final time of flight
%         - structure GNC, defining the characteristics of the algorithm

% Output: - array S, the rendezvous relative trajectory
%         - array u, the needed control vector 
%         - scalar tf, the final TOF 
%         - array lissajous_constants, indicating the evolution in time of
%           the Lissajous curve parameters

% New versions: implement DeCasteljau's algorithm

function [S, u, tf, lissajous_constants] = LSB_guidance(mu, L, gamma, s0, method, tf, GNC)    
    % Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end

    % Branch the different algorithms 
    switch (method)
       case 'Prescribed shape-based'
            [S, u, tf, lissajous_constants] = prescribed_lissajous(mu, L, gamma, s0, tf);
        case 'Numerical shape-based'
             [S, u, tf, lissajous_constants] = num_lissajous(mu, L, gamma, s0, GNC);
        case 'Dynamics shape-based'
            [S, u, tf, lissajous_constants] = dyn_lissajous(mu, L, gamma, s0, tf, GNC);
        case 'Minimum energy'
            [S, u, tf, lissajous_constants] = minimum_energy(mu, L, gamma, s0, tf, GNC);
        otherwise 
            error('No valid algorithm has been selected')
    end
end

%% Auxiliary functions
% Employ a 3rd order Bézier curve to impose the needed lissajous trajectory
function [S, u, tf, lissajous_constants] = prescribed_lissajous(mu, L, gamma, s0, tf)
    % Constants of the model
    cn = legendre_coefficients(mu, L, gamma, 2);                % Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 % Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 % In-plane frequency
    wv  = sqrt(c2);                                             % Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 % Constraint on the planar amplitude

    % Time span 
    tspan = 0:1e-3:tf;                                          % Dimensional time span 
    tau = tspan/tf;                                             % Non-dimensional time span

    % Initial constants
    phi0 = atan2(s0(2), -kap*s0(1));                            % In-plane initial phase
    lissajous_constants(7,:) = phi0*ones(1,length(tspan));      % In-plane initial phase 
    lissajous_constants(8,:) = (pi/2)*ones(1,length(tspan));    % Out-of-plane initial phase 

    % Amplitude computation 
    initial = -s0(1)/cos(phi0);                                 % Initial conditions on the in-plane amplitude
    initial([2 4]) = [s0(3); s0(6)];                            % Initial conditions on the amplitudes
    initial(3) = (-s0(4)+initial(1)*wp*sin(phi0))/cos(phi0);    % Initial conditions on the in-plane amplitude velocity
    final = zeros(4,1);                                         % Rendezvous sufficient/necessary conditions
    
    % Compute the Bézier curve polynomial basis 
    B = cell(2,1);                                                     
    for i = 1:2
        B{i} = [bernstein_basis(3,tau); bernstein_derivative(3,tau,1); bernstein_derivative(3,tau,2)];
    end

    initial(3:4) = tf*initial(3:4);                             % Non-dimensional initial amplitude velocities
    final(3:4) = tf*final(3:4);                                 % Non-dimensional final amplitude velocities
               
    % Control points
    P(:,1) = initial(1:2);
    P(:,2) = initial(1:2)+initial(3:4)/3;
    for i = 1:2
        P(i,3) = final(i)-final(length(final)/2+i)/3;
        P(i,4) = final(i);
    end

    % Bézier curve
    for i = 1:size(P,1)
        % State vector fitting
        for j = 1:3
            lissajous_constants(i+size(P,1)*(j-1),:) = P(i,1:4)*B{i}(1+4*(j-1):4*j,:);
        end
    end

    % Dimensionalisation 
    lissajous_constants(3:4,:) = lissajous_constants(3:4,:)/tf;
    lissajous_constants(5:6,:) = lissajous_constants(5:6,:)/tf^2;

    % Final TOF 
    tf = tspan(end);

    % Rendezvous phase-space trajectory
    Ax = lissajous_constants(1,:);                                           % In-plane amplitude
    Az = lissajous_constants(2,:);                                           % Out-of-plane amplitude
    phi = lissajous_constants(7,:);                                          % In-plane initial phase
    psi = lissajous_constants(8,:);                                          % In-plane initial phase
    dAx = lissajous_constants(3,:);                                          % First derivative of the in-plane amplitude
    dAz = lissajous_constants(4,:);                                          % First derivative of the out-of-plane amplitude
    ddAx = lissajous_constants(5,:);                                         % Second derivative of the in-plane amplitude
    ddAz = lissajous_constants(6,:);                                         % Second derivative of the out-of-plane amplitude
    
    S(1,:) = -Ax.*cos(wp*tspan+phi);                                         % X relative coordinate
    S(2,:) = kap*Ax.*sin(wp*tspan+phi);                                      % Y relative coordinate
    S(3,:) = Az.*sin(wv*tspan+psi);                                          % Z relative coordinate
    S(4,:) = -dAx.*cos(wp*tspan+phi)+wp*Ax.*sin(wp*tspan+phi);               % Vx relative velocity
    S(5,:) = kap*dAx.*sin(wp*tspan+phi)+kap*wp*Ax.*cos(wp*tspan+phi);        % Vy relative velocity
    S(6,:) = dAz.*sin(wv*tspan+psi)+wv*Az.*cos(wv*tspan+psi);                % Vz relative velocity
 
    % Compute the final control law 
    u = [-ddAx.*cos(wp*tspan+phi)+2*dAx*(wp-kap).*sin(wp*tspan+phi); ddAx*kap.*sin(wp*tspan+phi)+2*dAx*(kap*wp-1).*cos(wp*tspan+phi); ddAz.*sin(wv*tspan+psi)+2*dAz*wv.*cos(wv*tspan+psi)];
end

% Numerical solution
function [S, u, tf, lissajous_constants] = num_lissajous(mu, L, gamma, s0, GNC)
    % Setup 
    setup = GNC.SBOPT.setup;                                    % Setup of the shape-based optimization

    % Constants of the model
    cn = legendre_coefficients(mu, L, gamma, 2);                % Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 % Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 % In-plane frequency
    wv  = sqrt(c2);                                             % Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 % Constraint on the planar amplitude

    % Initial conditions
    phi0 = atan2(s0(2), -kap*s0(1));                            % Constant initial in-plane phase
    psi0 = pi/2;                                                % Constant initial out-of-plane phase

    % Amplitude computation 
    initial = -s0(1)/cos(phi0);                                 % Initial conditions on the in-plane amplitude
    initial([2 4]) = [s0(3); s0(6)];                            % Initial conditions on the out-of-plane amplitude and velocity
    initial(3) = (-s0(4)+initial(1)*wp*sin(phi0))/cos(phi0);    % Initial conditions on the in-plane amplitude velocity

    % Numerical optimization
    [x, ~, u, tf, ~, tau, ~, ~] = mfls_optimization(wp, wv, kap, phi0, psi0, initial, GNC.Tmax, setup); 

    % Rendezvous relative trajectory in phase space
    Ax = x(1,:);                                                             % In-plane amplitude
    Az = x(2,:);                                                             % Out-of-plane amplitude
    phi = repmat(phi0,1,size(x,2));                                          % In-plane initial phase
    psi = repmat(psi0,1,size(x,2));                                          % In-plane initial phase
    dAx = x(3,:);                                                            % First derivative of the in-plane amplitude
    dAz = x(4,:);                                                            % First derivative of the out-of-plane amplitude
    tspan = tf*tau;                                                          % Integration span
    S(1,:) = -Ax.*cos(wp*tspan+phi);                                         % X relative coordinate
    S(2,:) = kap*Ax.*sin(wp*tspan+phi);                                      % Y relative coordinate
    S(3,:) = Az.*sin(wv*tspan+psi);                                          % Z relative coordinate
    S(4,:) = -dAx.*cos(wp*tspan+phi)+wp*Ax.*sin(wp*tspan+phi);               % Vx relative velocity
    S(5,:) = kap*dAx.*sin(wp*tspan+phi)+kap*wp*Ax.*cos(wp*tspan+phi);        % Vy relative velocity
    S(6,:) = dAz.*sin(wv*tspan+psi)+wv*Az.*cos(wv*tspan+psi);                % Vz relative velocity

    lissajous_constants = x;
    lissajous_constants(7:8,:) = [phi; psi];
end

% Rendezvous as a dynamical system
function [S, u, tf, lissajous_constants] = dyn_lissajous(mu, L, gamma, s0, tf, GNC)
    % Constants of the model
    cn = legendre_coefficients(mu, L, gamma, 2);                % Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 % Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 % In-plane frequency
    wv  = sqrt(c2);                                             % Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 % Constraint on the planar amplitude

    % Initial conditions
    phi0 = atan2(s0(2), -kap*s0(1));                            % Constant in-plane phase
    psi0 = pi/2;                                                % Constant out-of-plane phase

    % Amplitude computation 
    initial = -s0(1)/cos(phi0);                                 % Initial conditions on the in-plane amplitude
    initial([2 4]) = [s0(3); s0(6)];                            % Initial conditions on the out-of-plane amplitude and velocity
    initial(3) = (-s0(4)+initial(1)*wp*sin(phi0))/cos(phi0);    % Initial conditions on the in-plane amplitude velocity
    final = zeros(4,1);                                         % Rendezvous sufficient/necessary conditions

    % Initial estimation of the TOF      
    tspan = 0:1e-3:tf; 

    % Setup of the integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', @(t,s)rendezvous(final, t, s)); 
    [t, x] = ode45(@(t,s)amplitude_dynamics(kap, wp, wv, phi0, psi0, t, s, GNC), tspan, initial);

    % Final TOF 
    tf = t(end); 

    % Rendezvous relative trajectory in phase-space
    Ax = x(:,1).';                                                           % In-plane amplitude
    Az = x(:,2).';                                                           % Out-of-plane amplitude
    phi = repmat(phi0,1,size(x,1));                                          % In-plane initial phase
    psi = repmat(psi0,1,size(x,1));                                          % In-plane initial phase
    dAx = x(:,3).';                                                          % First derivative of the in-plane amplitude
    dAz = x(:,4).';                                                          % First derivative of the out-of-plane amplitude

    S(1,:) = -Ax.*cos(wp*tspan+phi);                                         % X relative coordinate
    S(2,:) = kap*Ax.*sin(wp*tspan+phi);                                      % Y relative coordinate
    S(3,:) = Az.*sin(wv*tspan+psi);                                          % Z relative coordinate
    S(4,:) = -dAx.*cos(wp*tspan+phi)+wp*Ax.*sin(wp*tspan+phi);               % Vx relative velocity
    S(5,:) = kap*dAx.*sin(wp*tspan+phi)+kap*wp*Ax.*cos(wp*tspan+phi);        % Vy relative velocity
    S(6,:) = dAz.*sin(wv*tspan+psi)+wv*Az.*cos(wv*tspan+psi);                % Vz relative velocity

    % Compute the control law
    u = zeros(3,size(x,1));         % Preallocation for speed
    for i = 1:size(x,1)
        [dS, u([1 3],i)] = amplitude_dynamics(kap, wp, wv, phi0, psi0, t(i), x(i,1:4).', GNC);
        x(i,5:6) = dS(3:4);
    end

    lissajous_constants = x.';
    lissajous_constants(7:8,:) = [phi; psi];

    % Compute the final control law
    u = [u(1,:); zeros(1,length(tspan)); u(3,:)];     

    switch (GNC.Algorithm)
        case 'SDRE'
        otherwise
            u([1 3],:) = (u([1 3],:)+[2*dAx*(wp-kap).*sin(wp*tspan+phi); 2*dAz*wv.*cos(wv*tspan+psi)]).*[-cos(wp*tspan+phi); sin(wv*tspan+psi)];
    end

    % Y axis control 
    u(2,:) = x(:,5)*kap.*sin(wp*tspan+phi)+2*dAx*(kap*wp-1).*cos(wp*tspan+phi);
end

% Minimum energy solution 
function [S, u, tf, lissajous_constants] = minimum_energy(mu, L, gamma, s0, tf, GNC)
    % Constants of the model
    cn = legendre_coefficients(mu, L, gamma, 2);                % Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 % Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 % In-plane frequency
    wv  = sqrt(c2);                                             % Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 % Constraint on the planar amplitude
    Tmax = GNC.Tmax;                                            % Maximum axial acceleration

    % Time span 
    tspan = 0:1e-3:tf;                                          % Dimensional time span 
    tau = tspan/tf;                                             % Non-dimensional time span

    % Initial constants
    phi0 = atan2(s0(2), -kap*s0(1));                            % Initial in-plane phase
    lissajous_constants(7,:) = phi0*ones(1,length(tspan));      % Initial in-plane phase 
    lissajous_constants(8,:) = (pi/2)*ones(1,length(tspan));    % Initial out-of-plane phase 

    % Amplitude computation 
    initial = zeros(4,1);                                       % Preallocation of initial conditions
    initial(1) = -s0(1)/cos(phi0);                              % Initial condition on the in-plane amplitude
    initial([2 4]) = [s0(3); s0(6)];                            % Initial condition on the out-of-plane amplitude and velocity
    initial(3) = (-s0(4)+initial(1)*wp*sin(phi0))/cos(phi0);    % Initial condition on the in-plane amplitude velocity
    
    % Compute the analytical solution
    initial(3:4) = tf*initial(3:4);                             % Velocity dimensionalizing
    c1 = -6*initial(1:2)-4*initial(3:4);                        % Constant of motion
    c2 = 12*initial(1:2)+6*initial(3:4);                        % Constant of motion
    
    lissajous_constants(1:2,:) = (1/6)*c2.*tau.^3+(1/2)*c1.*tau.^2+initial(3:4)*tau+initial(1:2);
    lissajous_constants(3:4,:) = (1/2)*c2.*tau.^2+c1.*tau+initial(3:4);
    lissajous_constants(5:6,:) = c1+c2*tau;
               
    % Dimensionalisation 
    lissajous_constants(3:4,:) = lissajous_constants(3:4,:)/tf;
    lissajous_constants(5:6,:) = lissajous_constants(5:6,:)/tf^2;

    Ax = lissajous_constants(1,:);                                           % In-plane amplitude
    Az = lissajous_constants(2,:);                                           % Out-of-plane amplitude
    phi = lissajous_constants(7,:);                                          % In-plane initial phase
    psi = lissajous_constants(8,:);                                          % In-plane initial phase
    dAx = lissajous_constants(3,:);                                          % First derivative of the in-plane amplitude
    dAz = lissajous_constants(4,:);                                          % First derivative of the out-of-plane amplitude
    ddAx = lissajous_constants(5,:);                                         % Second derivative of the in-plane amplitude
    ddAz = lissajous_constants(6,:);                                         % Second derivative of the out-of-plane amplitude
    
    % Rendezvous relative trajectory in phase-space
    S(1,:) = -Ax.*cos(wp*tspan+phi);                                         % X relative coordinate
    S(2,:) = kap*Ax.*sin(wp*tspan+phi);                                      % Y relative coordinate
    S(3,:) = Az.*sin(wv*tspan+psi);                                          % Z relative coordinate
    S(4,:) = -dAx.*cos(wp*tspan+phi)+wp*Ax.*sin(wp*tspan+phi);               % Vx relative velocity
    S(5,:) = kap*dAx.*sin(wp*tspan+phi)+kap*wp*Ax.*cos(wp*tspan+phi);        % Vy relative velocity
    S(6,:) = dAz.*sin(wv*tspan+psi)+wv*Az.*cos(wv*tspan+psi);                % Vz relative velocity
 
    % Compute the final control law 
    u = Tmax*[-(ddAx+2*dAx*(wp-kap).*sin(wp*tspan+phi)).*cos(wp*tspan+phi); (ddAz+2*dAz*wv.*cos(wv*tspan+psi)).*sin(wv*tspan+psi)];
    u = [u(1,:); zeros(1,length(tspan)); u(2,:)];   
    u(2,:) = ddAx*kap.*sin(wp*tspan+phi)+2*dAx*(kap*wp-1).*cos(wp*tspan+phi);
end

% First order form of the amplitude dynamics 
function [dS, u] = amplitude_dynamics(kap, wp, wv, phi0, psi0, t, s, GNC)    
    % Constants
    Tmax = GNC.Tmax;               % Maximum thrust
    u = zeros(2,1);                % Preallocation of the control vector

    % Compute the control law 
    switch (GNC.Algorithm)
        case 'SDRE'
            % SDRE definition
            Q = GNC.LQR.StateMatrix;                            % State error weight matrix
            R = GNC.LQR.ControlMatrix;                          % Control effort weight matrix

            % Compute the amplitude vector field 
            Alpha = [2*(wp-kap)*tan(wp*t+phi0) 0; 0 -2*wv/tan(wv*t+psi0)];   % Acceleration dynamics
            Theta = [-1/cos(wp*t+phi0) 0; 0 1/sin(wv*t+psi0)];               % Control input acceleration matrixç

            % In-plane motion
            A = [0 1; 0 Alpha(1,1)];                                         % Constant state dynamics 
            B = [0; Theta(1,1)];                                             % Control input matrix
            K = lqr(A,B,Q,R);                                                % LQR matrix 
            u(1) = -K*s([1 3]);                                              % Control law

            % Out-of-plane motion
            A = [0 1; 0 Alpha(2,2)];                                         % Constant state dynamics 
            B = [0; Theta(2,2)];                                             % Control input matrix
            K = lqr(A,B,Q,R);                                                % LQR matrix 
            u(2) = -K*s([2 4]);                                              % Control law

           % Global dynamics 
           A = [zeros(2) eye(2); zeros(2) Alpha]; 
           B = [zeros(2); Theta];

        case 'Integrator SDRE'
            % SDRE definition
            Q = GNC.LQR.StateMatrix;                            % State error weight matrix
            R = GNC.LQR.ControlMatrix;                          % Control effort weight matrix

            % In-plane motion
            A = [0 1; 0 0];                                     % Double integrator dynamics
            B = [zeros(1); 1];                                  % Control input matrix
            K = lqr(A,B,Q,R);                                   % LQR matrix 
            u(1) = -K*s([1 3]);                                 % Control law

            % Out-of-plane motion
            B = [zeros(1); 1];                                  % Control input matrix
            K = lqr(A,B,Q,R);                                   % LQR matrix 
            u(2) = -K*s([2 4]);                                 % Control law

            % Compute the amplitude vector field 
            A = [zeros(2) eye(2); zeros(2,4)];                  % Constant state dynamics 
            B = [zeros(2); eye(2)];                             % Global control input matrix

        case 'Minimum time'
            e = Tmax*s(1)+(1/2)*abs(s(3))*s(3);                 % Error to the first switching curve
            u(1) = -Tmax*tanh(5e6*e);              
            e = Tmax*s(2)+(1/2)*abs(s(4))*s(4);                 % Error to the second switching curve
            u(2) = -Tmax*tanh(5e6*e);

            % Compute the amplitude vector field 
            A = [zeros(2) eye(2); zeros(2,4)];                  % Constant state dynamics 
            B = [zeros(2); eye(2)];                             % Global control input matrixç

        otherwise 
            error('No valid control algorithm was selected');
    end

    % Control input saturation 
    if (norm(u) > Tmax)
        u = Tmax*u/norm(u);
    end

    if (isnan(norm(u)))
        u = zeros(size(u));
    end

    % Vectorfield
    dS = A*s+B*u;
end

% Rendezvous event 
function [Pos, isterminal, dir] = rendezvous(final, ~, s)
    Pos = norm(final-s);        % Distance to the prescribed final state
    isterminal = 1;             % Halt the integration
    dir = 0;                    % Crossing direction of the event
end