%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/01/25
% File: Continuator.m 
% Issue: 0 
% Validated: 

%% Continuator %% 
% This class implements a generic continuator object %

classdef Continuator
    properties
        DiffCorrector;          % Differential corrector to be used 
        Config;                 % Configuration of the continuator
        ContinuationFunction;   % Constraint function to comply with
    end
    
    methods
        function [obj] = Continuator()
        end

        % Configuration 
        function [obj] = Configure(obj, ConfigOptions)
            obj.Config = ConfigOptions;
        end

        % Setter functions 
        function [obj] = set.ContinuationFunction(obj, myConstraints)
            if ( isa(myConstraints, 'function_handle') )
                obj.ContinuationFunction = myConstraints;
            else
                error('The constraint function shall be a function handle. Aborting...');
            end
        end

        % Pre-defined correctors 
        [FinalObject, Stats] = SingleParameterContinuation(obj, InitialGuess);  % Single parameter continuation
        [FinalObject, Stats] = PseudoArcContinuation(obj, InitialGuess);        % Pseudo-arc length continuation
    end
end

