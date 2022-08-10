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

    tol = 1e-5; 
    switch (manifold)
        case 'Center'
            alpha = zeros(6,size(C,2));
            S = cylindrical2cartesian(C(1:6,:),true);
            STM = stm_computation(mu, tf, St, n, P, sampling_distribution, basis, tau);
            for i = 1:size(alpha,2)
                aux = reshape(STM(:,i), [6 6]);
                [V,~] = eig(aux);
                alpha(:,i) = V^(-1)*S(1:6,i);
            end
            c = [c abs(alpha(3:6,:))-tol];
        case 'Unstable'
            alpha = zeros(6,size(C,2));
            S = cylindrical2cartesian(C(1:6,:),true);
            STM = stm_computation(mu, tf, St, n, P, sampling_distribution, basis, tau);
            for i = 1:size(alpha,2)
                aux = reshape(STM(:,i), [6 6]);
                [V,~] = eig(aux);
                alpha(:,i) = V^(-1)*S(1:6,i);
            end
            c = [c abs(alpha(1,end))-tol];
        case 'Stable'
            alpha = zeros(6,size(C,2));
            S = cylindrical2cartesian(C(1:6,:),true);
            STM = stm_computation(mu, tf, St, n, P, sampling_distribution, basis, tau);
            for i = 1:size(alpha,2)
                aux = reshape(STM(:,i), [6 6]);
                [V,~] = eig(aux);
                alpha(:,i) = V^(-1)*S(1:6,i);
            end
            c = [c abs(alpha(2,end))-tol];
        otherwise
    end
end