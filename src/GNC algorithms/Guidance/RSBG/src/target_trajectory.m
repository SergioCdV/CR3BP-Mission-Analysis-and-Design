%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 19/05/22

%% Target trajectory %%
% Function to compute the target's orbit trajectory in the non-dimensional
% time

% Inputs: - vector tspan, the integration time span 
%         - scalar T, the target's orbit period
%         - array Cp, the target's orbit Chebyshev polynomial weights

% Outputs: - array C, the final target's periodic evolution

function [C] = target_trajectory(tspan, T, Cp)
    % Preallocation of the Chebyshev polynomials
    P = zeros(size(Cp,2),length(tspan));
    
    % Periodicity check    
    if (tspan(end) > T)
        k = floor(tspan(end)/T);
        for i = 1:k
            index = find(tspan >= T);
            tspan(index(1):end) = tspan(index(1):end)-T;
        end
    end

    tspan = 2*tspan/T-1;
    for i = 1:length(tspan)
        P(:,i) = chebyshev('first', size(Cp,2), tspan(i));
    end   

    % Evaluate the target trajectory 
    C = Cp*P;
end