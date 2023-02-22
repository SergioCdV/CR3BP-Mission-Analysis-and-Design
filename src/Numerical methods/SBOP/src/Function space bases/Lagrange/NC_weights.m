%% Project:  
% Sergio Cuevas del Valle
% Date: 08/11/22
% File: NC_nodes
% Issue: 0 
% Validated: 

%% Newton-Cotes nodes %%
% This scripts provides the function to compute the Newton-Cotes nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the degree of the Lagrange polynomial of interest
%         - vector tau, the equally spaced domain in which the function is
%           evaluated

% Output: - vector y, the NC wweights

function [w] = NC_weights(N)
    % Prellocation 
    w = zeros(N+1,1);   
    dumb = linspace(0,N,1000);

    % Compute the Newton-Cotes weights 
    for i = 1:N+1
        % Compute the functional 
        func = ones(1,length(dumb));
        for j = 0:N 
            if ((i-1) ~= j)
                func = func .* (dumb-j);
            end
        end

        w(i) = (-1)^(N-(i-1))/(factorial(N-(i-1))*factorial((i-1)))*trapz(dumb,func);
    end
end