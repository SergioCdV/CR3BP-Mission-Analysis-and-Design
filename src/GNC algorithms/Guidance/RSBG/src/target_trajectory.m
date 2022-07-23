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
    % Construct the polynomial basis
%     if (tspan(1) ~= 0)
%         tspan = (tspan/tspan(end) + 1) * tspan(end)/2;
%     end
    
    % Periodicity check
    if (tspan(end) > T)
        k = floor(tspan(end)/T);
        if (tspan(1) ~= 0)
            index = zeros(1,k);
            for i = 1:k
                t = 2*i*T-tspan(end);
                vec = find(tspan >= t);
                index(i) = vec(1);
            end  

            tspan = tspan/tspan(end);
            for i = 1:k-1
                tspan(index(i):index(i+1)) = 2*tspan(index(i):index(i+1))/(tspan(index(i+1))-tspan(index(i)))-1;
            end
        else
            for i = 1:k
                index = find(tspan >= T);
                tspan(index(1):end) = tspan(index(1):end)-T;
            end
            tspan = 2*tspan/T-1;                          % Mapping into the T * [-1,1] domain
        end
    end
    
    P = CH_basis('first', size(Cp,2)-1, tspan);           % Chebyshev polynomials
    C = Cp*P;                                             % Evaluate the target trajectory
end