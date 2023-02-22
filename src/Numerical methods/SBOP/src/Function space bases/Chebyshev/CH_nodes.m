%% Project:  
% Sergio Cuevas del Valle
% Date: 07/05/20
% File: CH_nodes
% Issue: 0 
% Validated: 

%% Chebyshev nodes %%
% This scripts provides the function to compute the Chebysehv nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the degree of the Chebyshev polynomial of interest

% Output: - vector y, the Chebyshev nodes of interest

function [y] = CH_nodes(N)
    % Chebyshev nodes 
    i = N:-1:0;
    y = cos(pi*i/N);
end