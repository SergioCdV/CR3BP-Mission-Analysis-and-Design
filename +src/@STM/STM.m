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
       Eigenvalues;         % Floquet multipliers of the STM
       Eigenvectors;        % Floque eigenvectors of the STM
       StabilityIndex;      % Henon stability index of the STM
   end

   methods 
       % Constructor 
       function [obj] = STM( myStateDim )
            % Basic properties of the STM object
            obj.StateDim = myStateDim;
            obj.StabilityIndex = 1; 
            obj.Phi = eye(obj.StateDim);

            obj.Eigenvalues =  zeros(obj.StateDim,1);
            obj.Eigenvectors = zeros(obj.StateDim,obj.StateDim);
       end

       % Setters 
       function [obj] = set.Phi( obj, myPhi )

           if ( size(myPhi,1) ~= obj.StateDim )
               PhiAux = reshape( myPhi, [], 1 );
               obj.Phi = reshape( PhiAux, obj.StateDim, [] );
           else
               obj.Phi = myPhi;
           end

           % Auxiliary results
           [obj.Eigenvalues, obj.Eigenvectors] = obj.EigenDecomposition( obj.Phi );
           obj.StabilityIndex = obj.HenonStabilityIndex( obj.Eigenvalues );
       end
   end

   methods (Static)
       [lambda, v] = EigenDecomposition(STM);   % Eigendecomposition of the STM
       [nu] = HenonStabilityIndex(lambda);      % Stability index of the STM
       [CG] = CauchyGreenTensor(STM);           % Compute the Cauchy-Green tensor from a STM       
   end
end