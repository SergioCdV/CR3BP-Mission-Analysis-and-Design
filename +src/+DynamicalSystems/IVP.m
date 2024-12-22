%% KSHLand %% 
% Author: Sergio Cuevas
% Date: 15/11/2024

%% Initial Value Problem class %% 
% This is a class defining initial value problems (IVP) %

classdef IVP
    properties
            System;         % Dynamical system of the problem 
            IC;             % Initial conditions of the problem 
            t0 = 0;         % Initial value of the independent variable
    end

    methods 
        % Basic constructor 
        function [obj] = IVP(mySystem, myICs, myt0)
            % Sanity checks 
            if ( ~isa(mySystem, "DynamicalSystems.HybridSystem") )
                error('The input dynamical system for the IVP is not supported... Aborting');
            else
                obj.System = mySystem;
            end

            if ( ~exist("myICs", "var") )
                obj.IC = zeros( obj.System.StateDim, 1 );
            elseif ( size(myICs, 1) ~= obj.System.StateDim )
                error('Initial conditions do not match the state space dimension... Aborting');
            else
                obj.IC = myICs;
            end

            if ( ~exist("myt0", "var") )
                obj.t0 = 0;
            else
                obj.t0 = myt0;
            end
        end
    end
end