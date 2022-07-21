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
%         - string method, the parameter distribution to be used

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, St, T, initial, final, n, x, B, basis, tau, method)
    % Extract the optimization variables
    P = reshape(x(1:end-2), [length(n), max(n)+1]);     % Control points
    tf = x(end-1);                                      % Final time of flight 
    N = floor(x(end));                                  % Optimal number of revolutions

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);

    % Evaluate the target periodic trajectory 
    switch (St.Field)
        case 'Relative'
            tspan = tf*tau;
            St.Trajectory = target_trajectory(tspan, St.Period, St.Cp);
    end

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Control input 
    u = acceleration_control(mu, St, C, tf, method);

    % Equalities 
    ceq = [];

    % Inequality (control authority)
    switch (method)
        case 'Regularized'
            c = [sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-(tf^2*T*r.^2.*ones(1,size(u,2)))]; 
        otherwise
            c = [sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-(tf^2*T*ones(1,size(u,2)))];
    end
end