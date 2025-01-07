%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/24
% File: HaloOrbit.m 
% Issue: 0 

%% Halo Orbit class %%
% This class defines the baisc object to describe Halo orbits in the
% CR3BP

classdef HaloOrbit < src.OrbitFamilies.LissajousOrbit
    
    properties
        Branch = src.OrbitFamilies.HaloBranches.Northern;
    end
    
    methods
        % Constructor
        function [obj] = HaloOrbit( mySystem, myPoint, myBranch )
            % Parent constructor
            obj@src.OrbitFamilies.LissajousOrbit( mySystem, myPoint );

            if ( isa(myBranch, "src.OrbitFamilies.HaloBranches") )
                obj.Branch = myBranch;
            end
        end

        [seed] = OrbitSeed(obj, Amp, theta, order, freq, kap);  % Compute an initial seed for the orbit
    end

    methods (Static)
        
    end
end

