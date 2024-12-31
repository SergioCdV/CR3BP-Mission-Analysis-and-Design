%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 23/12/24
% File: DynamicsCR3BP.m 
% Issue: 0 
% Validated: 

%% CR3BP Dynamics %%
% This function contains the vector field of the CR3BP system. It accounts for a infinitesimal mass
% moving in the normalized, non dimensional synodic frame define by the two primaries, which
% are assumed to be in the same plane and in circular orbits. It also
% contains the integration of the first variational equations of the flow

% Inputs: 

% Outputs: - vector ds, the differential vector field of the system

% New versions: 

function [ds] = DynamicsCR3BP(obj, t, j, s, u, params)
    % Re-shaping of s 
    s = reshape(s, obj.StateDim, []);
    
    % Equations of motion of the CR3BP
    model = params{1};

    switch (model)
        % Deterministic models
        case "Newton"    
            ds = src.Systems.CR3BPSystem.NewtonEquationsCR3BP(t, j, s, u, params{2});       
            
        case "Encke"
            % Pre-allocation 
            ds = zeros(size(s));
            
            for i = 1:size(s,2)
                rel_state = s(1:obj.StateDim,i) - [obj.LP.r; zeros(3,size(obj.LP.r,2))];
                dist = sqrt( dot(rel_state(1:3,:), rel_state(1:3,:), 1) );

                rel_state = [obj.LP.r(:, dist == min(dist)); rel_state(:, dist == min(dist))];
    
                ds(:,i) = src.Systems.CoCR3BPSystem.EnckeEquationsCoCR3BP(t, j, rel_state, u(:,i), params{2});    
            end

        case "OrderN"
            ds = src.Systems.CoCR3BPSystem.NOrderEquationsCoCR3BP(t, j, s, u, params{2});    

        otherwise
            error('No valid CR3BP dynamics model has been selected. Aborting...');
    end
end