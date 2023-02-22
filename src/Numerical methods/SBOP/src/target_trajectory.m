%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 19/05/22

%% Target trajectory %%
% Function to compute the target's orbit trajectory in the non-dimensional
% time

% Inputs: - string sampling_distribution, indicating the distribution of the sampling grid
%         - vector tspan, the integration time span
%         - scalar T, the target's orbit period
%         - array Cp, the target's orbit Chebyshev polynomial weights

% Outputs: - array C, the final target's periodic evolution

function [C] = target_trajectory(t, T, Cp)
    theta = mod((2*pi/T)*t,2*pi);                         % Anomaly evolution
    tspan = 2*theta/(2*pi)-1;                             % Anomaly evolution
    P = CH_basis('first', size(Cp,2)-1, tspan);           % Chebyshev polynomials
    C = Cp*P;                                             % Evaluate the target trajectory
end