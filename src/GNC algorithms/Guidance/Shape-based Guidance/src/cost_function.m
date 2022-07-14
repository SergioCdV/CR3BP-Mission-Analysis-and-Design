%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Cost function %%
% Function to compute the cost function to be minimized

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector tau, the vector of collocation points
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used 
%         - string basis, the polynomial basis to be used
%         - string method, the parameter distribution to be used

% Outputs: - scalar r, the cost index to be optimized

function [r] = cost_function(mu, initial, final, n, tau, x, B, basis, method)
    % Minimize the control input
    P = reshape(x(1:end-2), [length(n), max(n)+1]);     % Control points
    tf = x(end-1);                                      % The final time of flight
    N = floor(x(end));                                  % The optimal number of revolutions
    N = 0;

    % Boundary conditions
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);

    % State evolution
    C = evaluate_state(P,B,n);

    % Control input
    u = acceleration_control(mu,C,tf,method);        

    % Control cost
    switch (method)
        case 'Regularized'
            r = sqrt(C(1,:).^2+C(3,:).^2);                   % Radial evolution
            a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);         % Non-dimensional acceleration
            a = a./r;                                        % Dimensional acceleration
        otherwise
            a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);         % Dimensional acceleration
    end
    
    % Cost function
    r = trapz(tau,a)/tf;                               
end