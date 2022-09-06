%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 21/05/21
% File: fcheb_derivative.m 
% Issue: 0 
% Validated: 

%% Derivative of a Chebyshev approximation %%
% This function allows the derivative of a function approximated by means of Chebyshev polynomials

% Inputs: - array C, the coefficients of the Chebyshev coefficients of the
%           function
%         - scalar b, the final value of the approximation domain
%         - scalar a, the initial value of the approximation domain

% Outpus: - vector dC, containing the evaluated derivative of the Chebyshev polynomials 

function [dC] = fcheb_derivative(C, b, a)
   % Sanity check on the order of the approximation
   n = size(C,2);                  
   
   % Preallocation of the Chebyshev polynomials and its derivatives
   dC = zeros(size(C));            % Derivative coefficients
   
   % Compute the derivative coefficients
   for i = 1:size(C,1)
       % Out-of-range cases
       dC(i,n) = 0;
       dC(i,n-1) = 2*(n-1)*C(i,n);
       
       % In-range cases
       for j = n-2:-1:1
          dC(i,j) = dC(i,j+2)+2*(j)*C(i,j+1);
       end
       
       dC(i,j) = dC(i,j)/2;        % Approximation term (C0 coefficient)
       
       % Renormalization
       for j = 1:size(C,2)
          dC(i,j) = (2/(b-a))*dC(i,j); 
       end
   end
end