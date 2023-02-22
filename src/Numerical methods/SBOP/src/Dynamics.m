%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Dynamics class %% 
% Class implementation of a dynamics function 

classdef Dynamics 
    % Fundamental definition of the cost function
    properties 
        path;           % Cost function file path 

        % Dynamics function
        ControlFunction;   
    end

    % Public methods
    methods 
        % Add cost function from txt file
        function [obj] = Dynamics()
            % Save the cost function 
            obj.ControlFunction = @(beta,t0,tf,x)myFunc(beta, t0, tf, s); 
        end
        
        % Evaluate the cost function
        function [u] = evaluate_control(obj, beta, t0, tf, s)
            % Evaluate the user defined cost function
            u = feval(myFunc, beta, t0, tf, s);
        end
    end
end