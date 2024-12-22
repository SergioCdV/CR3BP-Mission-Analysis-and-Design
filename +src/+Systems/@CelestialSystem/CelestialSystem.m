%% CR3BP Mission Design Library %% 
% Sergio Cuevas del Valle
% Date: 14/06/24
% File: CelestialSystem.m 
% Issue: 0 
% Validated:

%% Celestial System %% 
% This class implements the definition of a two body system, under whose
% influential gravity a third spacecraft moves (in the ballistic regime)

classdef (Abstract) CelestialSystem < src.DynamicalSystems.ContinuousSystem

    properties
        % Gravitational constants
        mu = 0;             % Reduced gravitational parameter of the system 
        M = zeros(1,2);     % Masses of the system [mass units]

        % Characteristic quantities
        Lc = 1;             % Characteristic length
        Fc = 1;             % Characteristic frequency 
        Tc;                 % Characteristic time
        Vc;                 % Characteristic velocity
        Ac;                 % Characteristic acceleration

        % System parameters 
        t0 = 0;             % Epoch of reference
        R = zeros(3,2);     % Position of the primaries in the synodic frame
    end

    methods
        % Constructor of the class
        % Inputs: - mu, scalar, reduced gravitational parameter of the
        %           system 
        %         - M, vector of 2x1, the gravitational masses of the
        %           system (first primary, secondary primary, in this order)
            
        % Output: - system object
        function [obj] = CelestialSystem(varargin)
           % Parent initialization 
           obj@src.DynamicalSystems.ContinuousSystem( 6, 0 );

           % Compute the reduced gravitational mass of the system
           if ( length(varargin) >= 2 )

               for i = 1:2:min(length(varargin), 2)
                   type = varargin{1 + 2 * (i-1)};
    
                   switch ( type )
                       case 'mu'
                           obj.mu = varargin{2};                % Mass of the system 
    
                       case 'M'
                           obj.M = varargin{2};                 % Masses of the primaries
                           obj.mu = obj.M(2) / sum(obj.M);      % Mass of the system
    
                       otherwise
                           error('Input parameters for the system are not supported. Aborting...');
                   end
               end

           elseif ( isempty(varargin) )
               error('Celestial system must have positive mass... Try again');

           else
               error('Input parameters for the system are not supported. Aborting...');

           end
        end

        % Function to initialize the rest of the parameters of the system 
        % Inputs:       - obj, the CelestialSystem object
        % Outputs:      - obj, the CelestialSystem object

        function [obj] = InitializeSystem(obj)
            % Compute the position of the primaries in Howell's synodic frame
            obj.R(1,:) = [-obj.mu 1-obj.mu];

            % Initialize characteristic quantities
            obj.Vc = 2*pi * obj.Lc / obj.Tc;            % Characteristic velocity
            obj.Ac = (2*pi)^2 * obj.Lc / obj.Tc^2;      % Characteristic acceleration
        end

        % Function to check if the masses of the primaries match its
        % parameters 
        % Inputs:       - obj, the CelestialSystem object
        % Outputs:      - check_flag, a boolean to output the result of the
        %                 test
        function [check_flag] = CheckSystem(obj)
            if ( abs( obj.M(2) / sum(obj.M) - obj.mu ) < 1E-5 )
                check_flag = true;
            else
                check_flag = false; 
                warning('The masses of the primaries do not correspond to the input system...');
            end
        end
    end
end