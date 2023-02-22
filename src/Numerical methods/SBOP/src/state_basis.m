%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 16/04/22

%% State basis %%
% Function to preallocate the state vector polynomial use of interest

% Inputs: - scalar n, the maximum polynomial order of the approximation 
%         - scalar L, the degree of the greatest derivative in the dynamics
%         - string basis, the polynomial basis to be used
%         - vector tau, the control parameter vector 

% Outputs: - cell array B, the state vector shape base basis

function [B, tau] = state_basis(n, L, basis, tau)
    % Preallocation of the Bernstein basis
    B = cell(length(n),1);              

    switch (basis)
        % Finite-horizon polynomial basis
        case 'Bernstein'
           for i = 1:length(n)
                B{i} = bernstein_basis(n(i),tau);
                for j = 1:L
                    B{i} = [B{i}; bernstein_derivative(n(i),tau,j)];
                end
           end
        case 'Orthogonal Bernstein'
            for i = 1:length(n)
                B{i} = OB_basis(n(i),tau);
                for j = 1:L
                    B{i} = [B{i}; OB_derivative(n(i),tau,j)];
                end
            end

        case 'Chebyshev'
            for i = 1:length(n)
                B{i} = CH_basis('first', n(i), tau);
                for j = 1:L
                    B{i} = [B{i}; CH_derivative('first',n(i),tau,j)];
                end
            end

        case 'Legendre'
             for i = 1:length(n)
                B{i} = LG_basis(n(i), tau);
                for j = 1:L
                    B{i} = [B{i}; LG_derivative(n(i),tau,j)];
                end
             end

        % Infinite-horizon polynomial basis
        case 'Hermite'
             for i = 1:length(n)
                B{i} = HT_basis(n(i), tau);
                for j = 1:L
                    B{i} = [B{i}; HT_derivative(n(i),tau,j)];
                end
             end
             
        case 'Laguerre'
             for i = 1:length(n)
                B{i} = LR_basis(n(i), tau);
                for j = 1:L
                    B{i} = [B{i}; LR_derivative(n(i),tau,j)];
                end
             end
        otherwise
            error('No valid functional polynomial basis has been selected');
    end
end