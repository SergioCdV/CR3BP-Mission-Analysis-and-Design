%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/01/25
% File: PeriodContinuator.m 
% Issue: 0 
% Validated: 

%% Period Continuator %% 
% This class implements a generic continuator along the period of a periodic object %

classdef PeriodContinuator < src.Continuation.Continuator
    properties
        Target;             % Target Jacobi Constant
    end
    
    methods
        function [obj] = PeriodContinuator( TargetPeriod )
            % Parent constructor 
            obj@src.Continuation.Continuator();

            % Save the target period 
            obj.Target = TargetPeriod;

            % Add the constraint function 
            obj.ContinuationFunction = @(InitialObject, ds)obj.TargetPeriodValue(obj.Target, InitialObject, ds);
        end
    end

    methods (Static)
        % Pre-defined constraint function
        function [e, NewObject] = TargetPeriodValue(TargetValue, InitialObject, ds)
            if ( isa(InitialObject, 'src.PeriodicOrbit') )
                % Compute the error 
                e = TargetValue - InitialObject.Period(1,1);

                % Update the orbit 
                NewObject = InitialObject;
                NewObject.Period = NewObject.Period + ds;
                
            else
                error('Continuation along the period of an object is only supported for periodic orbits. Aborting...');
            end
        end
    end
end

