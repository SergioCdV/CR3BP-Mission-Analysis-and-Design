%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/12/24
% File: Monodromy.m 
% Issue: 0 
% Validated: 

%% Monodromy Matrix %%
% This class definition provides the definition of several properties and
% methods associated to a Monodromy matrix

classdef Monodromy < src.STM
   properties
       Period;                  % Period of the monodromy matrix
       FloquetMultipliers;      % Floquet multipliers of the monodromy matrix 
       FloquetVectors;          % Floquet directions of the monodromy matrix
   end

   methods 
       % Constructor 
       function [obj] = Monodromy( myStateDim, myPeriod )
            % Parent constructor 
            obj@src.STM( myStateDim );

            % Local properties
            obj.Period = myPeriod;
       end
   end

   methods (Static)
       [f, Vf] = FloquetAnalysis(lambda, v, T);   % Transformation between Floquet multipliers and eigenvectors/eigenvalues
      [lambda] = LyapunovExponent();              % Compute the Lyapunov exponent of the STM 
   end
end