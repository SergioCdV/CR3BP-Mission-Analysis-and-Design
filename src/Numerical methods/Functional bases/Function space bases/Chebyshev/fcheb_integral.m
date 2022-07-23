%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/10/21
% File: fcheb_integral.m 
% Issue: 0 
% Validated: 

%% Integral of a Chebyshev approximation %%
% This function allows the integral of a function approximated by means of Chebyshev polynomials.

% Inputs: - array C, the coefficients of the Chebyshev coefficients of the
%           function
%         - scalar b, the final value of the approximation domain
%         - scalar a, the initial value of the approximation domain

% Outpus: - vector dC, containing the evaluated integral of the Chebyshev polynomials 

function [iC] = fcheb_integral(C, b, a)
   %Sanity check on the order of the approximation
   n = size(C,2);                  %Sanity check on the dimension
   
   %Preallocation of the Chebyshev polynomials and its integrals
   iC = zeros(size(C));            %Preallocation of the integrated coefficients
   fac = 1;                        %Sign coefficient 
   K = 0;                          %Constant of integration
   
   %Compute the derivative coefficients
   for i = 1:size(C,1)       
       %In-range cases
       for j = 2:n-1
          iC(i,j) = (C(i,j-1)-C(i,j+1))/(2*j-2);        %Coefficient
          K = K + fac*iC(i,j);
          fac = -fac; 
       end
       iC(i,n) = C(i,n-1)/(2*n-2); 
       K = K + fac*iC(i,end); 
       iC(i,1) = K; 
       iC(i,1) = iC(i,1);                             %Approximation term (C0 coefficient)

       %Renormalization
       for j = 1:size(C,2)
          iC(i,j) = (2/(b-a))*iC(i,j); 
       end
   end
end