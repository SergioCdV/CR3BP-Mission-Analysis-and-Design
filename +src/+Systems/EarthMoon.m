%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/12/24
% File: EarthMoon.m 
% Issue: 0 
% Validated:

%% CR3BP System %% 
% This class implements the definition of the Earth-Moon system

classdef EarthMoon < src.Systems.CR3BPSystem

    methods
        % Constructor of the class
            
        % Output: - system object
        function [obj] = EarthMoon(varargin)

           % System characteristics
           mu = 0.0121505856;                     % Mass parameter for the Earth-Moon system

           % Constructor of the super class
           obj@src.Systems.CR3BPSystem( 'mu', mu );

           % Characteristic quantities
           obj.M = [5.972e24 7.349e22];           % Masses of the system [kg]
           obj.Lc = 384399e3;                     % Mean distance from the Earth to the Moon [m]
           obj.Tc = 2.361e6;                      % Mean period of the Moon around the Earth [s]
           obj.Fc = 1 / obj.Tc;                   % Mean frequency of the Moon around the Earth [Hz]

           obj.CheckSystem(); 
           obj = obj.InitializeSystem();
        end
    end
end