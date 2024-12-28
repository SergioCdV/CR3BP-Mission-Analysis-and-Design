%% KSHLand %% 
% Author: Sergio Cuevas
% Date: 16/11/2024

%% Continuous System class %% 
% This is the basic, abstract class implementing a continuous system

classdef ContinuousSystem < src.DynamicalSystems.HybridSystem 
    properties
    end

    % Basic (abstract) methods
    methods
        % Constructor 
        function [obj] = ContinuousSystem(myStateDim, myParamsDim)
            % Sanity checks 
            if ( myStateDim <= 0 )
                warning('Input state dimension shall be strictly greater than 0...')
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
            FlowSet = @(t, j, x, u, params)( 1 );
            JumpSet = @(t, j, x, u, params)( 0 );
            JumpMap = @(t, j, x, u, params)( obj.NoJump(t, j, x, u) );

            obj.FlowSet = FlowSet;
            obj.JumpSet = JumpSet; 
            obj.Jump = JumpMap;
        end

        % Default jump mode
        function [t, x] = NoJump(obj, t, j, x, u)
            % Do nothing
        end
    end
end