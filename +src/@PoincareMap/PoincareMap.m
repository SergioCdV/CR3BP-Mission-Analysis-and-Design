%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/12/24
% File: PoincareMap.m 
% Issue: 0 
% Validated: 

%% Poincaré Map %%
% This class defines an abstract Poincare map object 

classdef PoincareMap
    properties
        NumberReturns;      % Integer number of allowed returns 
        SurfaceSection;     % Definition of the surface of section
    end
    
    methods
        % Constructor
        function [obj] = PoincareMap( myNumberReturns, mySoS )
            % Basic properties of the object
            obj.SurfaceSection = mySoS;
            obj.NumberReturns = myNumberReturns;
        end

        % Setters 
        function [obj] = set.SurfaceSection(obj, mySection)
            if isa(mySection, "function_handle")
                obj.SurfaceSection = mySection;
            else
                error('The input Surface of Section need be a function handle to continue...');
            end
        end

        function [obj] = set.NumberReturns(obj, myReturnNumber)
            if ( mod(myReturnNumber, 1) ~= 0 )
                warning('The input number of returns shall be integer...');
                myReturnNumber = floor(myReturnNumber);
            end

            obj.NumberReturns = myReturnNumber;
        end
        
        % Compute the Poincaré map
        [map, OrbitCollection] = Compute(obj, Solver, ICs);
    end
end

