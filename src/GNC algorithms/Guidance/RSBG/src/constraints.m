%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - scalar T, the maximum acceleration allowed for the spacecraft
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - vector tau, the vector of collocation points
%         - string sampling_distribution, the sampling distribution to be used
%         - string dynamics, the dynamics vectorfield parametrization to be
%           used
%         - string manifold, the manifold vector component to be nullified

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(cost, mu, St, T, initial, final, n, x, B, basis, tau, sampling_distribution, dynamics, manifold)
    % Extract the optimization variables
    P = reshape(x(1:end-2), [length(n), max(n)+1]);     % Control points
    tf = x(end-1);                                      % Final time of flight 
    N = floor(x(end));                                  % Optimal number of revolutions

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);

    % Evaluate the target periodic trajectory 
    switch (St.Field)
        case 'Relative'
            St.Trajectory = target_trajectory(sampling_distribution, tf, tau, St.Period, St.Cp);
    end

    % Trajectory evolution
    C = evaluate_state(P,B,n);
    S = cylindrical2cartesian(C(1:6,:), true);

    % Control input 
    [u, dv] = acceleration_control(mu, St, C, tf, dynamics);

    % Equalities
    switch (cost)
        case 'Minimum power'
            ceq = trapz(tau, dot(dv,u,1));
        otherwise
            ceq = [];
    end

    % Inequality (control authority)
    switch (dynamics)
        case 'Sundman'
            c = [sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-(tf^2*T*r.^2.*ones(1,size(u,2)))]; 
        case 'Euler'
            c = [sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-(tf^2*T*ones(1,size(u,2)))];
        otherwise
            error('No valid dynamics formulation was selected');
    end

    % Stable manifold departure constraint 
    Tr = manifold;                                                  % Relative period
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
    tspan = 0:1e-3:Tr;                                              % Integration time span
    s0 = [[target_trajectory(sampling_distribution, tf, 0, St.Period, St.Cp); target_trajectory(sampling_distribution, tf, 0, St.Period, St.Cv)]+C(1:6,1); reshape(eye(6), [36 1])];           % Relative initial conditions
    [~, Sr] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);

    STM = reshape(Sr(2,7:end), [6 6]);                                      % State transition matrix 
    [V, ~] = eig(STM);                                                      % Eigenspectrum of the STM 
 
    us = V(:,1)/norm(V(:,1));                                               % Stable manifold unit vector
    eps = 1e-4;                                                             % Stable manifold displacement
    s0(1:6) = s0(1:6)+eps*us;                                               % Stable relative initial conditions

    tspan = tf*tau;
    switch (sampling_distribution)
        case 'Chebyshev'
            tspan = 2*tspan;
            tspan = (tspan+2*tf)/2;
        case 'Legendre'
            tspan = 2*tspan;
            tspan = (tspan+2*tf)/2;
        otherwise
    end

    tspan = tspan(tspan < tf/10);
    tspan = flip(tspan);                                                                    % Reverse time span
    
    St = [target_trajectory(sampling_distribution, tf, flip(tspan/tf), St.Period, St.Cp); target_trajectory(sampling_distribution, tf, flip(tspan/tf), St.Period, St.Cv)];
   
    S = S(1:6,1:size(St,2))+St;                                                             % Relative trajecotry

    [~, Sr] = ode113(@(t,s)cr3bp_equations(mu, true, true, t, s), tspan, s0, options);      % Stable fiber

    alpha = sqrt(dot(Sr(:,1:3).'-S(1:3,:), Sr(:,1:3).'-S(1:3,:), 1));                       % Distance to the stable fiber
    c = [c alpha(1:end-1)-alpha(2:end)];                                                    % Monotinically departing away from the original stable fiber
    
end