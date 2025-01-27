%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/01/25
% File: JacobiContinuator.m 
% Issue: 0 
% Validated: 

%% Jacobi Continuator %% 
% This class implements a generic continuator along the Jacobi Constant %

classdef JacobiContinuator < src.Continuation.Continuator
    properties
        Target;             % Target Jacobi Constant
    end
    
    methods
        function [obj] = JacobiContinuator( TargetJC )
            % Parent constructor 
            obj@src.Continuation.Continuator();

            % Save the target Jacobi Constant 
            obj.Target = TargetJC;

            % Add the constraint function 
            obj.ContinuationFunction = @(InitialObject, ds)obj.TargetJacobiValue(obj.Target, InitialObject, ds);
        end
    end

    methods (Static)
        % Pre-defined constraint function
        function [e, NewObject] = TargetJacobiValue(TargetValue, InitialObject, ds)
            if ( isa(InitialObject, 'src.PeriodicOrbit') )
                % Update the orbit 
                NewObject = InitialObject;

                if ( isa(InitialObject, "src.OrbitFamilies.LyapunovOrbit") )
                    NewObject.State(1,1) = NewObject.State(1,1) + ds;

                elseif ( isa(InitialObject, "src.OrbitFamilies.HaloOrbit") )
                    NewObject.State(1,1) = NewObject.State(1,1) + ds;               % Moving along x to avoid the XZ bifurcation
                end

                % Compute the error 
                e = TargetValue - NewObject.EnergyFunction(1,1);
                
            else
                error('Continuation along the Jacobi Constant is only supported for periodic orbits. Aborting...');
            end
        end
    end
end

