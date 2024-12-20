%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 17/06/24
% File: CR3BPSystem.m 
% Issue: 0 
% Validated:

%% CR3BP System %% 
% This class implements the definition of a two body system, under whose
% influential gravity a third spacecraft moves (in the ballistic regime)

classdef CR3BPSystem < src.Systems.CelestialSystem

    properties
        % Additional properties of the system 
        LP;         % Libration points of the system 
    end

    methods
        % Constructor of the class
        % Inputs: - mu, scalar, reduced gravitational parameter of the
        %           system 
        %         - M, vector of 2x1, the gravitational masses of the
        %           system (first primary, secondary primary, in this order)
            
        % Output: - system object
        function [obj] = CR3BPSystem(varargin)
           % Constructor of the super class
           obj@src.Systems.CelestialSystem( varargin{:} );

           % Compute the libration points of the system 
           obj.LP = src.Systems.CR3BPSystem.LibrationPoints(obj.mu, obj.R);
        end
    end

    methods (Static)
        [Lp] = LibrationPoints(mu, R);              % Function to compute the libration points of the system
        [theta] = TimeLaw(T, t0, t);                % Time law describing the motion of the system
        [T] = Synodic2Inertial(theta);              % Homogeneous matrix (4x4) to transform from the inertial to the synodic reference frame
        [T] = Kepler2Synodic(idx, theta);           % Homogeneous matrix (4x4) to transform from the synodic barycentric to the synodic reference frame centered at one of the primaries
        
        [J, H] = JacobiConstant(mu, s);                         % Jacobi constant of the system 
        [U] = AugmentedPotential(mu, r);                        % Augmented potential of the system 
        [r] = ZeroVelocitySurface(mu, C, display_flag);         % Compute the ZVS associated to a given energy C 
        [r] = ZeroVelocityCurve(mu, C, display_flag);           % Compute the ZVC associated to a given energy C 
        [s] = ComplementaryZeroSurface(mu, C, display_flag);    % Compute the ZVC associated to a given energy C
    end

    methods (Access = private)
    end
end