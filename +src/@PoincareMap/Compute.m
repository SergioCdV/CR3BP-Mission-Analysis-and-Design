%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/12/24
% File: PoincareMap.m 
% Issue: 0 
% Validated: 

%% Compute Poincaré Map %%
% This script defines the function to compute a given Poincaré map

% Inputs: - HybridSolver Solver, the solver to integrate the map 
%         - cell array ICs, an array of Orbit objects defining the initial conditions to be
%           integrated (a blob of phase space)

% Output: - array map, the collection of phase space points on the Map
%           Surface of Section
%         - cell array OrbitCollection, an array of Orbit objects
%           containing the integrated trajectories

function [map, OrbitCollection] = Compute(obj, Solver, ICs, tspan)
    % Constants 
    N = size(ICs, 2);                           % Number of orbits to be computed

    % Pre-allocation 
    map = [];
    OrbitCollection = {};
    k = 1; 

    % Main loop 
    for i = 1:N
        try
            % Create the auxiliary hybrid system whose orbits are to be solved for 
            AuxSystem = ICs{i}.System;              % Copy of the original dynamical system 
           
            Solver.int_options = odeset(Solver.int_options, 'Events', @(t, j, s, u, params)( obj.SurfaceSection(t, j, s, u, params) ));
    
            % Define the IVP 
            IVP = src.DynamicalSystems.IVP( AuxSystem, ICs{i}.State(:,1), ICs{i}.t(1) );
    
            % Sanity checks 
            if ( ~exist("tspan", "var") )
                tspan = [IVP.t0 2*pi 1E-2];
            end
    
            ParticularSolver = Solver.configure( IVP );
            jspan = [0 obj.NumberReturns];
    
            % Solve the IVP
            [t, j, y, stats] = ParticularSolver.solve( tspan, jspan );
    
            % Check if the initial conditions or the second step are on the SoS
            for J = 1:2
                u = AuxSystem.ExogenousInput(t(J), j(J), y(:,J), AuxSystem.params);
                init_map = abs( obj.SurfaceSection(t(J), j(J), y(:,J), u, AuxSystem.params) ) < 1E-9;
                
                if (init_map)
                    map = [map y(:,J)];
                    break;
                end
            end
    
            % Check if the end conditions are on the SoS 
            if ( ~isempty(stats.EventValue) )
                event_states = [stats.EventValue{:}];
                map = [map event_states(3:end,:)];
            end
            
            trajectory{1} = t; 
            trajectory{2} = y;
            ICs{i}.State = trajectory;
            OrbitCollection{k} = ICs{i};
            k = k + 1;
        end
    end
end