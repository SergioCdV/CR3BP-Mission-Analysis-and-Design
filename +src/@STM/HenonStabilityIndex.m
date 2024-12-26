%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/12/24
% File: HenonStabilityIndex.m 
% Issue: 0 

%% Henon Stability Index %%
% This script contains the function to compute the Henon stability index
% associated with a dynamical solution of the problem

% Inputs: - array lambda [nxN], the eigenvalues are to be analyzed

% Output: - vector s [3xN], containing information about close bifurcations around the solution associated with the STM

% New versions: use symplecticity to correct computational error

function [s] = HenonStabilityIndex(lambda)
    % Sanity checks and constants 
    N = size(lambda, 2);                 % State dimension 
    
    % Pre-allocation 
    s = zeros(3, N);
    
    for i = 1:N
        lambda_aux = lambda(:,i);
        idx = imag(lambda_aux) == 0;

        % Neutrally stable Floquet multipliers
        if ( any(~idx) )
            s(3,i) = 0.5 * max( lambda_aux(~idx) + 1 ./ lambda_aux(~idx) );
        end

        % Unstable Floquet multipliers
        uns_idx = idx & ( abs(lambda_aux) > 1 );
        if ( any(uns_idx) )
            s(1,i) = 0.5 * max( lambda_aux(uns_idx) + 1 ./ lambda_aux(uns_idx) );
        end

        % Stable Floquet multipliers
        sta_idx = idx & ( abs(lambda_aux) < 1 );
        if ( any(sta_idx) )
            s(2,i) = 0.5 * max( lambda_aux(sta_idx) + 1 ./ lambda_aux(sta_idx) );   
        end
    end
end
