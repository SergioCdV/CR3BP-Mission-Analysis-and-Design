%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/01/25
% File: SunEarthMoon.m 
% Issue: 0 
% Validated:

%% CR3BP System %% 
% This class implements the definition of the Sun-Earth-Moon system

classdef SunEarthMoon < src.Systems.CR3BPSystem

    methods
        % Constructor of the class
            
        % Output: - system object
        function [obj] = SunEarthMoon(varargin)

           % System characteristics
           mu = 3.036E-6;                     % Mass parameter for the Sun-Earth-Moon system

           % Constructor of the super class
           obj@src.Systems.CR3BPSystem( 'mu', mu );

           % Characteristic quantities
           obj.M = [5.972e24 7.349e22];           % Masses of the system [kg]
           obj.Lc = 1.496E8;                      % Mean distance from the Earth-Moon to the Sun [m]
           obj.Tc = 3.147E7;                      % Mean period of the Earth-Moon around the Sun [s]
           obj.Fc = 1 / obj.Tc;                   % Mean frequency of the Earth-Moon around the Sun [Hz]

           obj.CheckSystem(); 
           obj = obj.InitializeSystem();
        end
    end
end