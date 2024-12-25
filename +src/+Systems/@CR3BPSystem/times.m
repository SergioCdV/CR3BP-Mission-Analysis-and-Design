%%  CR3BP Library %% 
% Author: Sergio Cuevas
% Date: 25/12/2024

%% Overloading of .* (times operator) %% 
% This scritp defines how the times operator (.*) is overloaded for
% CR3BPSystem, with the objective of defining a composition operator
% between the two inputs

function [OutSystem] = times(System1, System2)
    if ( isa(System2, "src.Systems.VariationalCR3BP") )
        % Pre-allocation
        OutSystem = System1;
        OutSystem.StateDim = System1.StateDim + System1.StateDim^2;        % Total state space dimension
    
        % Dynamics
        OutSystem.Dynamics = @(t, j, s, u, params)[OutSystem.Dynamics(t, j, s(1:System1.StateDim,:), u, params); System2.Dynamics(t, j, s, u, [System1.StateDim; params{2}])];
    else
        warning('The two input systems cannot be concatenated...')
    end
end