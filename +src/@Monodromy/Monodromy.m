%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/12/24
% File: Monodromy.m 
% Issue: 0 
% Validated: 

%% Monodromy Matrix %%
% This class definition provides the definition of several properties and
% methods associated to a Monodromy matrix

classdef Monodromy < src.Systems.STM
   properties
       Period;                  % Period of the monodromy matrix
       FloquetMultipliers;      % Floquet multipliers of the monodromy matrix 
       FloquetVectors;          % Floquet directions of the monodromy matrix
   end

   methods 
       % Constructor 
       function [obj] = Monodromy( mySTM, myPeriod )
            % Parent constructor 
            obj@src.Systems.STM( mySTM.StateDim );

            % Local properties
            obj.Period = myPeriod;
            obj.Phi = mySTM.Phi;
       end
   end

   methods (Static)
       [f, Vf] = FloquetAnalysis(lambda, v, T);   % Transformation between Floquet multipliers and eigenvectors/eigenvalues
      [lambda] = LyapunovExponent();              % Compute the Lyapunov exponent of the STM 
   end
end