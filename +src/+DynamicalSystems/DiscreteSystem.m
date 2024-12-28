%% KSHLand %% 
% Author: Sergio Cuevas
% Date: 16/11/2024

%% Continuous System class %% 
% This is the basic, abstract class implementing a continuous system

classdef DiscreteSystem < src.DynamicalSystems.HybridSystem 
    properties
    end

    % Basic (abstract) methods
    methods
        % Constructor 
        function [obj] = DiscreteSystem(myStateDim, myParamsDim)
            % Sanity checks
            if ( myStateDim <= 0 )
                warning('Input state dimension shall be strictly greater than 0...')
                return
            else
                StateDim = myStateDim;
            end

            if ( ~exist('myParamsDim', 'var') )
                ParamsDim = 0;
            else
                ParamsDim = myParamsDim;
            end

            % Parent class 
            obj@src.DynamicalSystems.HybridSystem(StateDim, ParamsDim);

            % Overall setting
            FlowSet = @(t, j, x, u, params)( 0 );
            JumpSet = @(t, j, x, u, params)( 1 );
            Dynamics = @(t, j, x, u, params)( x );

            obj.FlowSet = FlowSet;
            obj.JumpSet = JumpSet; 
            obj.Dynamics = Dynamics;
        end
    end
end