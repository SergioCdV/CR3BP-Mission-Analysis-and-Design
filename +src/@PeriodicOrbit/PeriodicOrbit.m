%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 29/12/24
% File: PeriodicOrbit.m 
% Issue: 0 

%% Periodic Orbit class %%
% This class defines the baisc object to describe a particular periodic solution of
% the CR3BP

classdef PeriodicOrbit < src.Orbit
    
    properties
        Period;         % Period of the periodic orbit
        Monodromy;      % Monodromy matrix of the periodic orbit 
    end
    
    methods
        % Constructor
        function [obj] = PeriodicOrbit( myStateDim, mySystem )
            obj@src.Orbit(myStateDim, mySystem);
        end
    end
end

