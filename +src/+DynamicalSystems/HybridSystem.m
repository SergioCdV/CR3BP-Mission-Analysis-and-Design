%% KSHLand %% 
% Author: Sergio Cuevas
% Date: 15/11/2024

%% Hybrid System class %% 
% This is the basic, abstract class implementing a hybrid dynamical system,
% defined by the flow F, jump J sets and their corresponding dynamics and
% jump map

classdef HybridSystem 
    properties
        StateDim;           % Dimension of the state space
        ParamsDim;          % Dimension of the parameter space
        ControlDim = 0;     % Dimension of the control space

        FlowSet;            % Handler for the flow set of the system 
        JumpSet;            % Handler for the jump set of the system 
        
        Dynamics;           % Dynamics in the flow set
        Jump;               % Set-value map in the jump set
        InputSignal;        % Handler for the controller

        PriorityRule = 1;   % Used to prioritize jumps or flows; default, prioritizes jumps

        params;
    end

    % Basic (abstract) methods
    methods
        % Constructor 
        function [obj] = HybridSystem(myStateDim, myParamsDim)
            % Sanity checks 
            if ( myStateDim <= 0 )
                warning('Input state dimension shall be strictly greater than 0...')
                return
            end

            obj.StateDim = myStateDim;

            if ( ~exist('myParamsDim', 'var') )
                obj.ParamsDim = 0;
            else
                obj.ParamsDim = myParamsDim;
            end

            % Default input signal
            obj.ControlDim = obj.StateDim;
            obj.InputSignal = @(t, j, x, params)( zeros(obj.ControlDim, 1) );
        end

        % Check if the state is in the flow set
        function [within_flag] = inFlowSet(obj, t, j, x, u, params)

            within_flag = obj.FlowSet( t, j, x, u, params );
            
        end

        % Check if the state is in the jump set
        function [within_flag] = inJumpSet(obj, t, j, x, u, params)

             within_flag = obj.JumpSet( t, j, x, u, params );

        end

        %% Setters 
        % Dynamics
        function [obj] = set.FlowSet(obj, myFlowSet)
            if ( ~isa(myFlowSet, 'function_handle') )
                error('The flow set needs to be a function handle... Aborting')
            else
                obj.FlowSet = myFlowSet;
            end
        end

        % Jump 
        function [obj] = set.JumpSet(obj, myJumpSet)
            if ( ~isa(myJumpSet, 'function_handle') )
                error('The jump set needs to be a function handle... Aborting')
            else
                obj.JumpSet = myJumpSet;
            end
        end

        % Dynamics
        function [obj] = set.Dynamics(obj, myDynamics)
            if ( ~isa(myDynamics, 'function_handle') )
                error('The dynamics map needs to be a function handle... Aborting')
            else
                obj.Dynamics = myDynamics;
            end
        end

        % Jump 
        function [obj] = set.Jump(obj, myJump)
            if ( ~isa(myJump, 'function_handle') )
                error('The jump map needs to be a function handle... Aborting')
            else
                obj.Jump = myJump;
            end
        end

        % Input signal 
        function [obj] = set.InputSignal(obj, myInputSignal)
            if ( ~isa(myInputSignal, 'function_handle') )
                error('The input signal map needs to be a function handle... Aborting')
            else
                obj.InputSignal = myInputSignal;
            end
        end

        % Overloaded operations 
        [OutSystem] = times(System2, System1);
    end

    methods(Static)
    end
end