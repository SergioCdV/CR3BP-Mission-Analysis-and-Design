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
%         - scalar tf, the final time of flight
%         - vector s0, initial conditions of both the target and the chaser
%         - scalar epsilon, the manifold displaement factor
%         - scalar tol, the differential corrector scheme tolerance for the
%           constrained maneuver

% Output: - array S, the rendezvous relative trajectory
%         - array u, the needed control vector 
%         - scalar tf, the final TOF 
%         - array lissajous_constants, indicating the evolution in time of
%           the Lissajous curve parameters

% New versions: 

%% Auxiliary functions 
function [s0, u, tf, lissajous_constants] = LSB_guidance(mu, s0, method, tf)    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end

    %Branch the different algorithms 
    switch (method)
        case 'Optimal shape-based'
        case 'Prescribed shape-based'
            [lissajous_constants, tspan] = prescribed_lissajous(mu, s0, tf);
        case 'SDRE'
        case 'Primer'
        otherwise 
            error('No valid algorithm has been selected')
    end

    % Final TOF 
    tf = tspan(end);
    
    % Compute the final relative trajectory 
    cn = legendre_coefficients(mu, L, gamma, 2);                %Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 %Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 %In-plane frequency
    wv  = sqrt(c2);                                             %Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 %Contraint on the planar amplitude

    Ax = lissajous_constants(1,:);        % In-plane amplitude
    Az = lissajous_constants(2,:);        % Out-of-plane amplitude
    phi = lissajous_constants(7,:);       % In-plane initial phase
    psi = lissajous_constants(8,:);       % In-plane initial phase
    dAx = lissajous_constants(3,:);       % First derivative of the in-plane amplitude
    dAz = lissajous_constants(4,:);       % First derivative of the out-of-plane amplitude
    ddAx = lissajous_constants(5,:);      % Second derivative of the in-plane amplitude
    ddAz = lissajous_constants(6,:);      % Second derivative of the out-of-plane amplitude
    
    S(1,:) = -Ax.*cos(wp*tspan+phi);                                         % X relative coordinate
    S(2,:) = kap*Ax.*sin(wp*tspan+phi);                                      % Y relative coordinate
    S(3,:) = Az.*sin(wv*tspan+psi);                                          % Z relative coordinate
    S(4,:) = -dAx.*cos(wp*tspan+phi)+wp*Ax.*sin(wp*tspan+phi);               % Vx relative velocity
    S(5,:) = kap*Ax.*sin(wp*tspan+phi)+kap*wp*Ax.*cos(wp*tspan+phi);         % Vy relative velocity
    S(6,:) = dAz.*sin(wv*tspan+psi)+wv*Az.*cos(wv*tspan+psi);                % Vz relative velocity

    % Compute the final control law 
    u = [ddAx.*cos(wp*tspan+phi); -kap*ddAx.*sin(wp*tspan+phi); ddAz.*sin(wv*tspan+psi)];
end

%% Auxiliary functions
% Employ a 4th order Bézier curve to impose the needed lissajous trajectory
function [lissajous_constants, tspan] = prescribed_lissajous(mu, L, gamma, s0, tf)
    % Constants of the model
    cn = legendre_coefficients(mu, L, gamma, 2);                %Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 %Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 %In-plane frequency
    wv  = sqrt(c2);                                             %Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 %Contraint on the planar amplitude

    % Time span 
    tspan = 0:1e-3:tf;                                          %Dimensional time span 
    tau = tspan/tf;                                             %Non-dimensional time span

    % Initial constants
    phi0 = atan2(s0(2), -kap*s0(1));
    lissajous_constants(7,:) = phi0*ones(1,length(tspan));      % In-plane initial phase 
    lissajous_constants(8,:) = pi/2*ones(1,length(tspan));      % Out-of-plane initial phase 

    % Amplitude computation 
    initial = [-s0(1)/cos(phi0); s0(3); (s0(4)-wp*sin(phi0))/cos(phi0); s0(6)];     % Initial conditions
    final = zeros(4,1);                                                             % Rendezvous sufficient conditions
    
    % Compute the Bézier curve 
    B = cell(2,1);                  % Polynomial basis           
    for i = 1:2
        B{i} = [bernstein_basis(3,tau); bernstein_derivative(3,tau,1); bernstein_derivative(3,tau,2)];
    end

    initial(3:4) = tf*initial(3:4);
    final(3:4) = tf*final(3:4);
               
    P(:,1) = initial(1:2);
    P(:,2) = initial(1:2)+initial(3:4)/3;
    for i = 1:2
        P(i,3) = final(i)-final(length(final)/2+i)/3;
        P(i,4) = final(i);
    end

    for i = 1:size(P,1)
        % State vector fitting
        for j = 1:3
            lissajous_constants(i+N*(j-1),:) = P(i,1:4)*B{i}(1+4*(j-1):4*j,:);
        end
    end
end