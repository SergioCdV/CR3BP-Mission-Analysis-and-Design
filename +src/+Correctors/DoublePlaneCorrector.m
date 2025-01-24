%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/01/25
% File: DoublePlaneCorrector.m 
% Issue: 0 
% Validated: 

%% Double Plane Differential Corrector %% 
% This class implements a differential corrector object for orbits doubly symmetric around the XZ and YZ planes %

classdef DoublePlaneCorrector < src.Correctors.DiffCorrector    
    methods
        % Basic constructor
        function [obj] = DoublePlaneCorrector()
            % Parent object
            obj@src.Correctors.DiffCorrector( 2, 2 )

            % Constraint function 
            obj.ComputeConFunction = @(Orbit, ds)obj.DoublyConstraint(Orbit, ds);
        end
    end

    methods (Static)
        % Pre-defined constraint function 
        function [e, M, Orbit] = DoublyConstraint(Orbit, ds)
            % Apply the computed correction 
            Orbit.State([1 6],1) = Orbit.State([1 6],1) + ds;

            % Ensure motion on the XY synodic plane
            Orbit.State(2,1) = 0;                           % Null Y coordinate 
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
            jspan = [0 2];
            max_event_cnt = 1;
            [t, j, y, ~] = Solver.solve( tspan, jspan, max_event_cnt );

            STM.Phi = y(Orbit.StateDim+1:end,end);          % Monodromy matrix
            Orbit.STM = STM;                                % STM of the system
            Orbit.Period = 4 * t(end);                      % New orbital period
            
            % Compute the error
            e = [y(4,end); y(6,end)];                       % Vx and Vz at the crossing must be 0
            
            % Compute the linear system matrix 
            Phi = STM.Phi;                                  % Build the monodromy matrix at T/2

            % Vector field at T
            u = zeros(Orbit.System.ControlDim,1);
            F = Orbit.System.Dynamics(t(end), j(end), y(1:Orbit.StateDim,end), u, Orbit.System.params);   

            % Final matrix
            M = [Phi(4,1) Phi(4,6); Phi(6,1) Phi(6,6)] - [F(4,1); F(6,1)] / y(5,end) * [Phi(2,1) Phi(2,6)];
        end
    end
end

