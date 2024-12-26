%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/12/24
% File: Orbit.m 
% Issue: 0 

%% Orbit class %%
% This class defines the baisc object to describe a particular solution of
% the CR3BP

classdef Orbit
    
    properties
        StateDim;           % Phase space dimension
        System;             % Dynamical system associated for which the orbit is a solution

        State;              % Phase space trajectory
        t;                  % Independent variable
        STM;                % State transition matrix along the trajectory

        EnergyFunction;     % Jacobi constant of the phase trajectory
    end
    
    methods
        % Constructor
        function [obj] = Orbit( myStateDim, mySystem )
            % Basic properties
            obj.StateDim = myStateDim;              % State dimension
            obj.System = mySystem;                  % Dynamical system 

            % Initial basic values 
            obj.t = 0;                              % Template initial independent variable for the trajectory
            obj.State = zeros(obj.StateDim,1);      % Template initial conditions
        end
        
        % Setter
        function [obj] = set.State(obj, myState)
            % Independent variable and sanity checks 
            if ( size(myState,2) == 2)
                obj.t = myState{1};
            else
                obj.t = [];
            end

            if isa(myState, "cell")
                myState = myState{end};
            end

            % Set the trajectory
            if ( size(myState,1) ~= obj.StateDim )
               StateAux = reshape( myState, [], 1 );
               obj.State = reshape( StateAux, obj.StateDim, [] );
            else
               obj.State = myState;
            end

            % Compute the Jacobi constant along the trajectory 
            if ( isa(obj.System, "src.Systems.CR3BPSystem") )
                [~, obj.EnergyFunction] = src.Systems.CR3BPSystem.JacobiConstant(obj.System.mu, obj.State);
            else
                obj.EnergyFunction = nan;
            end
        end

        function [obj] = set.STM(obj, mySTM)
            if ( isa(mySTM, "src.Systems.STM") )
                obj.STM = mySTM;
            else
                warning('The input State Transition Matrix shall be an STM object...');
                obj.STM = [];
            end
        end
    end
end

