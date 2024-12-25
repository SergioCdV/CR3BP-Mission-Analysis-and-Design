%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/12/24
% File: HenonStabilityIndex.m 
% Issue: 0 

%% Henon Stability Index %%
% This script contains the function to compute the Henon stability index
% associated with a dynamical solution of the problem

% Inputs: - array STM [nxnN], whose eigenvalues are to be analyzed

% Output: - vector nu [1xN], containing information about close bifurcations around the solution associated with the STM

% New versions: use symplecticity to correct computational error

function [s] = HenonStabilityIndex(STM)
    % Sanity checks and constants 
    n = size(STM,1);                    % State dimension 

    if ( mod(size(STM,2), n) == 0 )
        N = size(STM,2) / n;            % Number of STM 
    
        % Pre-allocation 
        s = zeros(1, N);
        
        for i = 1:N
            idx = 1 + n * (i-1) : n * i;

            % Compute the eigenspectrum of the STM 
            [V, D] = eig( STM(:,idx) );
                       
            % Henon stability indices
            if (flag)
                s(1) = (1/2) * ( eig(1,1)+eig(6,6) );    % Sum of the reciprocal pair
                s(2) = (1/2) * ( eig(2,2)+eig(3,3) );    % Sum of the neutrally stable pair
                s(3) = (1/2) * ( eig(4,4)+eig(5,5) );    % Sum of the remaining pair
            else
                s(i) = 0;
            end
        end
    else
        warning('The second dimension of the STM does not match the expected value...');
        s = 0;
    end
end
