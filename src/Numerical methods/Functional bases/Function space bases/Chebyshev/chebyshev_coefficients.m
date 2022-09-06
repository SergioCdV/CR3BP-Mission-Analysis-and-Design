%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 21/05/21
% File: chebyshev_coefficients.m 
% Issue: 0 
% Validated: 

%% Chebyshev polynomials approximation %%
% This function allows to compute the coefficients of a Chebyshev expansion
% to approximate a given one-dimensional function

% Inputs: - data domain x
%         - data y, the vector to regress or approximate 
%         - scalar order, the order of the approximation

% Outpus: - vector cn, containing the coefficients of the approximation 
%         - matrix W, the matrix of weights to scale the approximation

function [cn, W] = chebyshev_coefficients(x, y, order)
    % Sanity check on the data dimensions
    if (size(y,1) == 1)
        y = y.';
    end
    
    % Constants 
    a = x(1);                                        % Initial point in the domain
    b = x(end);                                      % Final point in the domain
    
    % Project the data on the interval [-1, 1]
    u = (2*x-(b+a))/(b-a);                           % New argument
    
    % Main computation 
    Tn = CH_basis('first', order, u).';              % Compute the Chebyshev polynomials
    
    % Sampling distribution characteristics
    W = diag(ones(1,length(y)));
    
    % Chebyshev coefficients by least-squares
    cn = (Tn.'*W*Tn)^(-1)*(Tn.'*W*y);
    cn = cn.';
end

