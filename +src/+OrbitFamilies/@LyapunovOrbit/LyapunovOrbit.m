%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/24
% File: LyapunovOrbit.m 
% Issue: 0 

%% Lyapunov Orbit class %%
% This class defines the baisc object to describe Lyapunov orbits in the
% CR3BP

classdef LyapunovOrbit < src.OrbitFamilies.LissajousOrbit
    
    properties
    end
    
    methods
        % Constructor
        function [obj] = LyapunovOrbit( mySystem, myPoint )
            % Parent constructor
            obj@src.OrbitFamilies.LissajousOrbit( mySystem, myPoint );

            % Nullify the out-of-plane component
            obj.OrbitFrequencies(2) = 0;
        end
    end
end

