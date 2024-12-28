%%  CR3BP Library %% 
% Author: Sergio Cuevas
% Date: 27/12/2024

%% Overloading of .* (times operator) %% 
% This scritp defines how the times operator (.*) is overloaded for
% CoCR3BPSystem, with the objective of defining a composition operator
% between the two inputs

function [OutSystem] = times(System1, System2)
    if ( isa(System2, "src.Systems.CR3BPSystem") )
        % Pre-allocation
        OutSystem = System1;
    
        OutSystem.PhaseSpaceDim = [System2.PhaseSpaceDim; System1.PhaseSpaceDim];
        OutSystem.OriginalStateDim = [System2.StateDim; System1.StateDim];
        OutSystem.StateDim = System1.StateDim + System2.StateDim;                       % Total state space dimension
        
        OutSystem.ParamsDim = 1; 
        OutSystem.params{2} = System2.params{2};

        OutSystem.VariationalProblem = [OutSystem.VariationalProblem; System2.VariationalProblem];
    
        % Dynamics
        OutSystem.Dynamics = @(t, j, s, u, params)[System2.Dynamics(t, j, s(1:System2.StateDim,:), System2.ExogenousInput(t, j, s(1:System2.StateDim,:), System2.params), System2.params); ...
                                                   OutSystem.DynamicsCoCR3BP(t, j, s, u, params)];

    elseif ( isa(System2, "src.Systems.VariationalCoCR3BP") )
        OutSystem = System1;
        
        if (System1.StateDim^2 ~= System2.StateDim)
            error('The dimension of the variational problem does not match that of the original system. Aborting...');
        else
            OutSystem.StateDim = System1.StateDim + System1.StateDim^2;        % Total state space dimension
        end
    
        % Dynamics
        OutSystem.Dynamics = @(t, j, s, u, params)[OutSystem.Dynamics(t, j, s(1:System1.StateDim,:), u, params); System2.Dynamics(t, j, s, u, [System1.StateDim; params{2}])];

        OutSystem.VariationalProblem(1) = true;

    else
        warning('The two input systems cannot be concatenated...')
    end
end