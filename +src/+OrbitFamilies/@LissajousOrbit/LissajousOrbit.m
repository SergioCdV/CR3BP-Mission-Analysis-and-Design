%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/24
% File: LissajousOrbit.m 
% Issue: 0 

%% Lissajous Orbit class %%
% This class defines the baisc object to describe Lissajous orbits in the
% CR3BP

classdef LissajousOrbit < src.PeriodicOrbit
    
    properties
        LibrationPoint;         % Index of the libration point
        Origin;                 % Location of the libration orbit
        OrbitFrequencies;       % Vertical and planar frequencies

        OrbitAmplitudes;        % Vertical and planar amplitudes
        OrbitPhases;            % Vertical and planar initial condition phases
    end

    properties (Hidden)
        kap;                    % Amplitude constraint
    end
    
    methods
        % Constructor
        function [obj] = LissajousOrbit( mySystem, myPoint )
            % Parent constructor
            obj@src.PeriodicOrbit(6, mySystem);

            % Additional properties 
            if ( myPoint > 0 && myPoint <= 5 )
                obj.LibrationPoint = myPoint;
                obj.Origin = mySystem.LP.r(:,obj.LibrationPoint);

                % Frequencies of the orbit
                [obj.OrbitFrequencies, obj.kap] = obj.LissajousFrequencies( mySystem.mu, obj.LibrationPoint, mySystem.LP.gamma(obj.LibrationPoint) );

                % Period of the orbit
                obj.Period = (2*pi) / obj.OrbitFrequencies(1);                                
            end
        end

        % Setters
        function [obj] = set.OrbitAmplitudes(obj, myAmp)
            if ( any(myAmp < 0) )
                error('The orbit amplitudes shall be stricitly positive. Aborting...')
            
            else
                obj.OrbitAmplitudes = myAmp;
            
            end
        end

        % Additional methods
        [seed] = OrbitSeed(obj, Amp, theta, freq, kap);      % Compute an initial seed for the orbit
        [phi, num_rev] = TimeLaw(obj, tspan, theta);         % Time law on the orbit
    end

    methods (Static)
        [w, kap] = LissajousFrequencies(mu, L, gamma);       % Compute the frequencies of the Lissajous orbit
    end
end

