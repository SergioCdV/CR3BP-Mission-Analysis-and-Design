%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/12/24
% File: CauchyGreenTensor.m 
% Issue: 0 

%% Cauchy-Green Tensor %%
% This script computes the Cauchy-Green tensor for a given STM

% Inputs: - array STM [nxnN], the STM whose CGT should be computed

% Output: - array CGT [nxnN], the evolution of the CGT to be computed

function [CGT] = CauchyGreenTensor(STM)
    % Sanity checks and constants 
    n = size(STM,1);                    % State dimension 

    if ( mod(size(STM,2), n) == 0 )
        N = size(STM,2) / n;            % Number of STM 
    
        % Pre-allocation 
        CGT = zeros(n, n * N);
        
        for i = 1:N
            idx = 1 + n * (i-1) : n * i;

            % Compute the CGT of the STM 
            CGT(:,idx) = STM(:,idx).' * STM(:,idx);
        end
    else
        warning('The second dimension of the STM does not match the expected value...');

        CGT = [];
    end
end