%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/09/22
% File: bernstein_basis.m 
% Issue: 0 
% Validated:

%% Bernstein basis %%
% Function for calculation the Bernstein polynomial basis of degree n

% Inputs:  - scalar n, the cardinal of the basis
%          - vector tvec, the control parameter t vector 

% Output: - array B, the Bernstein polynomial basis of dimensions    
%           n x length(t)

% Iterative portion of De Casteljau's algorithm
function [B] = bernstein_basis(n, tvec)    
    % Find number of steps (time increments)
    steps = length(tvec);
    
    % Initialize variable for n-order curve
    B = zeros(n+1,steps); 
    
    % Calculation of bezier curve, recursively
    for i = 0:n
        % Indexing variable 
        j = i+1; 
        
        % binomial coefficient (n i)
        bin = factorial(n)/(factorial(i)*factorial(n-i));
        
        % Bernstein basis polynomial
        B(j,:) = bin.*tvec.^i.*(1-tvec).^(n-i);
    end
end