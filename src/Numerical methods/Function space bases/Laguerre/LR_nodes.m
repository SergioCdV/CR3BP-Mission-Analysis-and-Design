%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/09/22
% File: LR_nodes.m 
% Issue: 0 
% Validated:

%% Laguerre nodes %%
% This scripts provides the function to compute the Laguerre nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the degree of the Legendre polynomial of interest
%         - scalar alpha, the generalized Laguerre polynomial coefficient

% Output: - vector x, the Laguerre nodes of interest

function [x] = LR_nodes(N, alpha)        
        % Map of the Legendre nodes to the Laguerre domain
        x = LG_nodes(N);
        x = (x+1)/2;
        x = tan(0.95*pi/2*x);
end