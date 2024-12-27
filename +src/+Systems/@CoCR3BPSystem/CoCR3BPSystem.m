%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: CoCR3BPSystem.m 
% Issue: 0 
% Validated:

%% Co-orbital CR3BP System %% 
% This class implements the definition of the co-orbital problem in a two body system, under whose
% influential gravity a third spacecraft moves (in the ballistic regime)

classdef CoCR3BPSystem < src.DynamicalSystems.ContinuousSystem

    properties
        % Additional properties of the system 
        ForceModel;         % Function handle for the perturbations model
        ControlInput;       % Control signal to the system 
    end

    properties (Access = private)
    end

    methods
        % Constructor of the class
        % Inputs: - mu, scalar, reduced gravitational parameter of the
        %           system 
        %         - M, vector of 2x1, the gravitational masses of the
        %           system (first primary, secondary primary, in this order)
            
        % Output: - system object
        function [obj] = CoCR3BPSystem(varargin)
           % Constructor of the super class
           myStateDim = 6;                                              % The state vector is 3 position + 3 velocity
           obj@src.DynamicalSystems.ContinuousSystem(myStateDim, 0);

           % Complete the system 
           obj.ControlDim = 3;      % The control dimension is 3 (acceleration)

           model = src.Systems.ModelsCoCR3BP.Newton;
           obj.params{1} = model;

           % Function handles 
           obj.ForceModel =   @(t, j, x, params)( zeros(obj.ControlDim, 1) );
           obj.ControlInput = @(t, j, x, params)( zeros(obj.ControlDim, 1) );
           obj.Dynamics =     @(t, j, s, u, params)obj.DynamicsCoCR3BP(t, j, s, u, params);
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
        [U] = AugmentedPotential(mu, r);                        % Augmented potential of the system 
        [J, H] = JacobiConstant(mu, s);                         % Jacobi constant of the system 

        [cn] = CoLegendreCoefficients(mu, r_t, order);              % Legendre coefficients of the co-orbital Hamiltonian

        [ds] = NewtonEquationsCoCR3BP(t, j, s, u, params);          % Newton's description of the co-orbital CR3BP dynamics
        [ds] = EnckeEquationsCoCR3BP(t, j, s, u, params);           % Encke's description of the co-orbital CR3BP dynamics
        [ds] = LinearEquationsCoCR3BP(t, j, s, u, params);          % Linear model of the co-orbital problem
        [ds] = SecondOrderEquationsCoCR3BP(t, j, s, u, params);     % Second order model of the co-orbital problem 
        [ds] = ThirdOrderEquationsCoCR3BP(t, j, s, u, params);      % Third order model of the co-orbital problem
        [ds] = LibrationEquationsCoCR3BP(t, j, s, u, params);       % Linear model of the co-orbital problem
        [ds] = RichardsonEquationsCoCR3BP(t, j, s, u, params);      % Linear model of the co-orbital problem around a collinear libration point
        [J] = JacobianCoCR3BP(mu, s);                               % Jacobian of the co-orbital dynamics vector field
    end

    methods (Access = private)
    end
end