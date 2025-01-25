%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/01/25
% File: PlanarCorrector.m 
% Issue: 0 
% Validated: 

%% Planar Differential Corrector %% 
% This class implements a differential corrector object for planar orbits %

classdef PlanarCorrector < src.Correctors.DiffCorrector    
    methods
        % Basic constructor
        function [obj] = PlanarCorrector()
            % Parent object
            obj@src.Correctors.DiffCorrector( 1, 1 )

            % Constraint function 
            obj.ComputeConFunction = @(Orbit, ds)obj.PlanarConstraint(Orbit, ds);
        end
    end

    methods (Static)
        % Pre-defined constraint function 
        function [e, M, Orbit] = PlanarConstraint(Orbit, ds)
            % Apply the computed correction 
            Orbit.State(5,1) =  Orbit.State(5,1) + ds;

            % Ensure motion on the XY synodic plane
            Orbit.State(2,1) = 0;                           % Null Y coordinate 
            Orbit.State(3,1) = 0;                           % Null Z coordinate
            Orbit.State(4,1) = 0;                           % Null Vx 
            Orbit.State(6,1) = 0;                           % Null Vz 
            
            % Augment initial conditions with the initial STM 
            STM = src.STM( Orbit.StateDim );
            STM.Phi = eye( Orbit.StateDim );

            % Create the complete variational system 
            VarSystem = src.Systems.VariationalCR3BP( Orbit.StateDim );
            CompleteSystem = Orbit.System .* VarSystem;

            % Integrator 
            options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14, 'Events', @(t, j, s, u, params)( s(2) ) );
            integrator = src.DynamicalSystems.HybridSolver( @ode113, options );

            % Initial Value Problem 
            s0 = [Orbit.State(:,1); reshape(STM.Phi, [], 1)];
            VarCR3BPIVP = src.DynamicalSystems.IVP( CompleteSystem, s0, Orbit.t );
            
            % Configuration 
            Solver = integrator.configure( VarCR3BPIVP );

            % Solve the system 
            tspan = [Orbit.t(1) 2*pi 0.01];
            [t, j, y, ~] = Solver.solve( tspan );

            STM.Phi = y(Orbit.StateDim+1:end,end);          % Monodromy matrix
            Orbit.STM = STM;                                % STM of the system
            Orbit.Period = 2 * t(end);                      % New orbital period
            
            % Compute the error
            e = y(4,end);                                   % Vx at the crossing must be 0
            
            % Compute the linear system matrix 
            Phi = STM.Phi;                                  % Build the monodromy matrix at T/2

            % Vector field at T/2
            u = zeros(Orbit.System.ControlDim,1);
            F = Orbit.System.Dynamics(t(end), j(end), y(1:Orbit.StateDim,end), u, Orbit.System.params);   

            % Final matrix
            M = Phi(4,5) - ( F(4,1) / y(5,end) ) * Phi(2,5);
        end
    end
end

