%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/04/22

%% Variational equations %%
% Function to compute the variational equations and the linear evolution of
% the STM

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the 9xm state vector
%         - scalar tf, the final time of flight
%         - string dynamics, the dynamics vectorfield formulation to be
%           used
%         - string sampling_distribution, the sampling distribution to be used

% Outputs: - vector dPhi, the linear variational equations 

function [dPhi] = var_equations(mu, tf, St, n, P, sampling_distribution, basis, t, s)
    % Reshape the STM
    Phi = reshape(s, [6 6]);      % STM 

    % Compute the chaser and target trajectories 
    rt = target_trajectory(sampling_distribution, tf, t, St.Period, St.Cp);        % Target trajectory
    B = state_basis(n, t, basis);                                                  % Approximation basis 
    C = evaluate_state(P, B, n);                                                   % Chaser trajectory
    rho = cylindrical2cartesian(C(1:6), true);                                     % Cartesian chaser trajectory 
    s = [rt; zeros(3,1); rho];                                                     % Complete relative state

    % Relative jacobian 
    A = rel_jacobian(mu, s);

    % Final variational equations 
    dPhi = A*Phi; 
    dPhi = reshape(dPhi, [36 1]);
end