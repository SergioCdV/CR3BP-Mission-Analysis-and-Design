%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 23/12/24
% File: STM.m 
% Issue: 0 
% Validated: 

%% State Transition Matrix %%
% This class definition provides the definition of several properties and
% methods associated to a State Transition matrix

classdef STM
   properties
       StateDim;            % Dimension of the STM
       Phi;                 % State transition matrix per se 
       StabilityIndex;      % Henon stability index of the STM
   end

   methods 
       % Constructor 
       function [obj] = STM( myStateDim )
            % Basic properties of the STM object
            obj.StateDim = myStateDim;
            obj.StabilityIndex = 1; 
            obj.Phi = eye(obj.StateDim);
       end

       % Setters 
       function [obj] = set.Phi( obj, myPhi )

           if ( size(myPhi,2) == 1 )
               obj.Phi = reshape( myPhi, obj.StateDim, [] );
           else
               obj.Phi = myPhi;
           end
           
           obj.StabilityIndex = obj.HenonStabilityIndex( obj.Phi );
       end
   end

   methods (Static)
       [nu] = HenonStabilityIndex(STM);     % Stability index of the STM
       [CG] = CauchyGreenTensor(STM);       % Compute the Cauchy-Green tensor from a STM
       [lambda] = LyapunovExponent();       % Compute the Lyapunov exponent of the STM 
       
   end
end