%% KSHLand %% 
% Author: Sergio Cuevas
% Date: 18/11/2024

%% Tsien Problem class %% 
% This is the basic, abstract class implementing the classical radial
% thrust problem, also known as Tsien's

classdef Tsien < DynamicalSystems.RegularizationKS.TwoBPKS
    methods
        % Constructor
        function [obj] = Tsien( myEps )
            % Parent constructor
            super_arguments{1} = 1;

            obj@DynamicalSystems.RegularizationKS.TwoBPKS( super_arguments{:} );

            % Parameters
            obj.ParamsDim = 2;
            obj.params(2) = myEps / 8;

            obj.InputSignal = @(t, j, s, params)obj.TsienInputSignal(t, j, s, params);
            obj.ControlDim = 3;
        end
    end 

    methods (Static)
        % Methods
        function [ap] = TsienInputSignal(t, j, s, params)
            % Cartesian state 
            X = DynamicalSystems.RegularizationKS.KS_mapping( s, false, "1", params(1) );

            % Oblateness perturbation 
            ap = params(2) * X(1:3,:) ./ sqrt( dot(X(1:3,:), X(1:3,:), 1) );
        end
    end
end