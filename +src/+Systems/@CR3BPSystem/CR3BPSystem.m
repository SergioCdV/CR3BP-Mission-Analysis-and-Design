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
        LP;                 % Libration points of the system 
        ForceModel;         % Function handle for the perturbations model
        ControlInput;       % Control signal to the system 
    end

    properties (Hidden)
        PhaseSpaceDim = 6;                  % Dimension of the phase space
        VariationalProblem = false;         % Boolean flag to indicate if the problem is augmented with the variational system 
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
           obj.R(:,1) = [-obj.mu; 0; 0];
           obj.R(:,2) = [1 - obj.mu; 0; 0];

           obj.LP = src.Systems.CR3BPSystem.LibrationPoints(obj.mu, obj.R);

           % Complete the system 
           obj.StateDim = 6;        % The statae vector is 3 position + 3 velocity
           obj.ControlDim = 3;      % The control dimension is 3 (acceleration)
           obj.ParamsDim = 7;       % The dynamics depend only on mu

           model = src.Systems.ModelsCR3BP.Newton;
           obj.params{1} = model;
           obj.params{2} = [obj.mu; reshape(obj.R, [], 1)];

           % Function handles 
           obj.ForceModel =   @(t, j, x, params)( zeros(obj.ControlDim, 1) );
           obj.ControlInput = @(t, j, x, params)( zeros(obj.ControlDim, 1) );
           obj.Dynamics =     @(t, j, s, u, params)obj.DynamicsCR3BP(t, j, s, u, params);
        end

        % Setters 
        function [obj] = set.ForceModel(obj, myForceModel)
            if ( ~isa(myForceModel, 'function_handle') )
                error('The force model function needs to be a function handle... Aborting')
            else
                obj.ForceModel = myForceModel;
                obj.ExogenousInput = @(t, j, s, params)( obj.ForceModel(t, j, s, params) + obj.ControlInput(t, j, s, params) );
            end
        end

        % Input signal 
        function [obj] = set.ControlInput(obj, myControlSignal)
            if ( ~isa(myControlSignal, 'function_handle') )
                error('The control signal map needs to be a function handle... Aborting')
            else
                obj.ControlInput = myControlSignal;
                obj.ExogenousInput = @(t, j, s, params)( obj.ForceModel(t, j, s, params) + obj.ControlInput(t, j, s, params) );
            end
        end
        
        % Propagation of the CR3BP dynamics
        [ds] = DynamicsCR3BP(obj, t, j, s, u, params);          % Vector field of the system
    end

    methods (Static)
        [theta] = TimeLaw(T, t0, t);                            % Time law describing the motion of the system
        [T] = Synodic2Inertial(theta, direction);               % Homogeneous matrix (4x4) to transform from the inertial to the synodic reference frame
        [T] = Kepler2Synodic(mu, idx, theta, direction);        % Homogeneous matrix (4x4) to transform from the synodic barycentric to the synodic reference frame centered at one of the primaries
        
        [Lp] = LibrationPoints(mu, R);                          % Function to compute the libration points of the system
    
        [U] = PotentialFunction(mu, r);                         % Potential function of the system 
        [U] = AugmentedPotential(mu, r);                        % Augmented potential of the system 
        [J, H] = JacobiConstant(mu, s);                         % Jacobi constant of the system 
        [r] = ZeroVelocitySurface(mu, C, display_flag);         % Compute the ZVS associated to a given energy C 
        [r] = ZeroVelocityCurve(mu, C, display_flag);           % Compute the ZVC associated to a given energy C 
        [s] = ComplementaryZeroSurface(mu, C, display_flag);    % Compute the ZVC associated to a given energy C

        [ds] = NewtonEquationsCR3BP(t, j, s, u, params);        % Newton's description of the CR3BP dynamics
        [J] = JacobianCR3BP(mu, s);                             % Jacobian of the absolute dynamics vector field

        [c] = LegendreCoefficients(mu, L, gamma, order);        % Legendre coefficients to expand the CR3BP Hamiltonian around the libration points
    end

    methods (Access = private)
    end
end