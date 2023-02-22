%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost class %% 
% Class implementation of a cost function 

classdef Cost 
    % Fundamental definition of the cost function
    properties 
        path;           % Cost function file path 

        cost;           % Value of the cost function
        units;          % Units of the cost function

        CostFunction;   % Cost function
    end

    % Public methods
    methods 
        % Add cost function from txt file
        function [obj] = Cost()
            % Save the cost function 
            obj.CostFunction = @(beta,t0,tf,x,u)CostFunction(beta, t0, tf, s, u); 
        end

        % Add units
        function [obj] = AddUnits(obj, myUnits)
            % Units of the cost function
            obj.units = myUnits;          
        end
        
        % Evaluate the cost function
        function [M, L] = evaluate_cost(obj, beta, t0, tf, s, u)
            % Evaluate the user defined cost function
            cost = feval(myFunc, beta, t0, tf, s, u);

            % Decompose the cost vector
            M = cost(1);    % Mayer penalty term
            L = cost(2);    % Lagrange integral term
        end
    end

    methods (Access = private)
        % Add cost function from txt file
        function [obj] = CostFromTxt(obj, myCostFile)
            % Read the txt file
            obj.path = strcat('Problem\',myCostFile);
            eqs = fileread(obj.path);
            eqs = strsplit(eqs, newline);
            eqs = strtrim(eqs);

            % Remove empty lines 
            Eqs = {};
            k = 1;
            for i = 1:size(eqs,2)
                if (~isempty(eqs{i}))
                    Eqs{k} = eqs{i};
                    k = k+1;
                end
            end

            myFunc = str2func(['@(beta,t0,tf,x,u) ' strcat(Eqs{1,:})]);

            % Save the cost function 
            obj.CostFunction = myFunc; 
        end
    end
end