%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/04/22
% File: CH_nodes.m 
% Issue: 0 
% Validated: 

%% Chebyshev nodes %%
% This scripts provides the function to compute the Chebysehv nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the degree of the Chebyshev polynomial of interest

% Output: - vector y, the Chebyshev nodes of interest

function [y] = CH_nodes(N)
    % Chebyshev nodes 
    i = N-1:-1:0;
    y = cos(pi*i/(N-1));
end