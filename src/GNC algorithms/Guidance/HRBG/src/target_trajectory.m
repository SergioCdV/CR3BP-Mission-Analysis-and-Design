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

function [C] = target_trajectory(sampling_distribution, tf, tau, T, Cp)
    % Construct the polynomial basis
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
    
    k = floor(tspan(end)/T);
    for i = 1:k
        index = find(tspan >= T);
        tspan(index(1):end) = tspan(index(1):end)-T;
    end

    tspan = 2*tspan/T-1;                                  % Mapping into the T * [-1,1] domain
    P = CH_basis('first', size(Cp,2)-1, tspan);           % Chebyshev polynomials
    C = Cp*P;                                             % Evaluate the target trajectory
end