%% Project:  
% Sergio Cuevas del Valle
% Date: 20/01/20
% File: HT_nodes
% Issue: 0 
% Validated: 

%% Hermite nodes %%
% This scripts provides the function to compute the Gauss-Hermite nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the degree of the Hermite polynomial of interest

% Output: - vector x, the Legendre nodes of interest

function [x] = HT_nodes(N)
    % Map of the Legendre nodes to the Laguerre domain
    x = LG_nodes(N);
    x = tan(0.95*pi/2*x);
end