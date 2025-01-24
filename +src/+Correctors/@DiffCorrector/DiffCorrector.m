%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/01/25
% File: DiffCorrector.m 
% Issue: 0 
% Validated: 

%% Differential Corrector %% 
% This class implements a generic differential corrector object %

classdef DiffCorrector
    properties
        Config;                 % Configuration options for the differential corrector
        ConDim;                 % Number of constraints of the problem
        DoFDim;                 % Number of degrees of freedom of the system 

        ComputeConFunction;     % Constraint function
    end
    
    methods
        % Basic constructor
        function [obj] = DiffCorrector( numVariables, numConstraints )
            obj.ConDim = numConstraints;
            obj.DoFDim = numVariables;
        end

        % Configuration 
        function [obj] = Configure(obj, ConfigOptions)
            obj.Config = ConfigOptions;
        end

        % Setter functions 
        function [obj] = set.ComputeConFunction(obj, myConstraints)
            if ( isa(myConstraints, 'function_handle') )
                obj.ComputeConFunction = myConstraints;
            else
                error('The constraint function shall be a function handle. Aborting...');
            end
        end
        
        % Pre-defined correctors 
        [FinalOrbit, Stats] = SingleShootSolve(obj, InitialGuess);          % Single shooting correction
        [FinalOrbit, Stats] = MultipleShootSolve(obj, InitialGuess);        % Multiple shooting correction
    end
end

