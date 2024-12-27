%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: DynamicsCoCR3BP.m 
% Issue: 0 
% Validated: 

%% Co-orbital CR3BP Dynamics %%
% This function contains the vector field of the co-orbital CR3BP system

% Inputs: 

% Outputs: - vector ds, the differential vector field of the system

% New versions: 

function [ds] = DynamicsCoCR3BP(obj, t, j, s, u, params)
    
    % Restrict the state to the right components 
    if ( obj.StateDim > 6 )
        if ( obj.StateDim == 12 )
            % Do nothing, the state is target + co-orbital state

        elseif ( obj.StateDim == 72 )
            % For the dynamics, s should be target + co-orbital state (no variational components)
            tgt_idx = 1:obj.StateDim;
            cor_idx = (obj.StateDim^2 + obj.StateDim + 1) : (obj.StateDim^2 + 2 * obj.StateDim);
            s = s([tgt_idx cor_idx],:);

        else
            % For the dynamics, s should be target + co-orbital state (no variational components)
            tgt_idx = obj.StateDim + 1 : 2 * obj.StateDim;
            cor_idx = 1:obj.StateDim : (obj.StateDim^2 + 2 * obj.StateDim);
            s = s([tgt_idx cor_idx],:);
        end
    end

    % Equations of motion of the CR3BP
    model = params{1};

    switch (model)
        % Deterministic models
        case "Newton"    
            ds = src.Systems.CoCR3BPSystem.NewtonEquationsCoCR3BP(t, j, s, u, params{2});   
            
        case "Encke"    
            ds = src.Systems.CoCR3BPSystem.EnckeEquationsCoCR3BP(t, j, s, u, params{2});    
          
        case "Linear"
            params{2} = [params{2}; 0];         % Do not consider 2nd order effects 

            ds = src.Systems.CoCR3BPSystem.LinearEquationsCoCR3BP(t, j, s, u, params{2});  

        case "Order2"
            params{2} = [params{2}; 0];         % Do not consider 3rd order effects 

            ds = src.Systems.CoCR3BPSystem.SecondOrderEquationsCoCR3BP(t, j, s, u, params{2});   

        case "Order3"
            params{2} = [params{2}; 0];         % Do not consider 4th order effects 

            ds = src.Systems.CoCR3BPSystem.ThirdOrderEquationsCoCR3BP(t, j, s, u, params{2}); 

        case "OrderN"
            params{2} = [params{2}; 0];         % Do not consider (N+1)-th order effects 

            error('The selected co-orbital dynamics model is not supported yet...');

        case "Libration"
            ds = src.Systems.CoCR3BPSystem.LibrationEquationsCoCR3BP(t, j, s, u, params{2}); 

        case "Richardson"
            ds = src.Systems.CoCR3BPSystem.RichardsonEquationsCoCR3BP(t, j, s, u, params{2});

        otherwise
            error('No valid CR3BP dynamics model has been selected. Aborting...');
    end
end