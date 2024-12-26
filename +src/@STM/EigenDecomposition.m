%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/12/24
% File: EigenDecomposition.m 
% Issue: 0 

%% Eigen-decomposition %%
% This script contains the function to compute the decomposition of the
% of the STM in eigenvalues and eigenvectors

% Inputs: - array STM [nxnN], whose eigenvalues are to be analyzed

% Output: - vector lambda [1xN], containing the eigenvalues of the STM
%         - array v [nxnN], containing the eigenvectors of the STM

function [lambda, v] = EigenDecomposition(STM)
    % Sanity checks and constants 
    n = size(STM,1);                    % State dimension 

    if ( mod(size(STM,2), n) == 0 )
        N = size(STM,2) / n;            % Number of STM 
    
        % Pre-allocation 
        lambda = zeros(n, N);
        v = zeros(n, n * N);
        
        for i = 1:N
            idx = 1 + n * (i-1) : n * i;

            % Compute the eigenspectrum of the STM 
            [V, D] = eig( STM(:,idx) );

            lambda(:,i) = diag( D );
            v(:,idx) = V; 
        end
    else
        warning('The second dimension of the STM does not match the expected value...');

        lambda = []; 
        v = [];
    end
end