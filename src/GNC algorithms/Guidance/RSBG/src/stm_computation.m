%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% State Transition Matrix Computation %% 
% Function to compute the STM of the problem

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - scalar tf, the final time of flight 
%         - structure St, defining the target trajectory in time
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - array P, the control points of the trajectory
%         - string sampling_distribution, the sampling distribution to be used
%         - string basis, the polynomial basis to be used
%         - vector tau, the vector of collocation points

% Outputs: - array STM, the evolution of the STM in time  


function [STM] = stm_computation(mu, tf, St, n, P, sampling_distribution, basis, tau)
        % Integration
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);                                            % Integration setup
        Phi0 = reshape(eye(6), [36 1]);                                                                   % Initial conditions

        [~, STM] = ode113(@(t,s)var_equations(mu, tf, St, n, P, sampling_distribution, basis, t, s), tau, Phi0, options);
        STM = STM.'; 
end