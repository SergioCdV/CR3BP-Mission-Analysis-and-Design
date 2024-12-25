%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/12/24
% File: CheckSymplecticity.m 
% Issue: 0 

%% Check Symplecticity %%
% This script contains the function to compute the Henon stability index
% associated with a dynamical solution of the problem

% Inputs: - array A [nxn], the matrix whose symplecticity is to be checked

% Output: - double error, the numerical error to check the symplecticity of
%           the input array

function [error] = CheckSymplecticity(A)
    % Constants 
    n = size(A,1);          % Dimension of the group

    if ( mod(n,2) == 0 )
        I = eye(n/2);
        O = zeros(n/2);
        J = [O I; -I O];    % Canonical symplectic matrix

        % Error 
        error = A.' * (J * A) - J;
        [~, S, ~] = svd(error, 'econ');
        error = max( max(S) );
    else
        error = Inf;        % The matrix cannot be symplectic by definition
    end    
end