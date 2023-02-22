%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Dynamics class %% 
% Class implementation of an Optimal Problem function 

classdef Problem 
    % Fundamental definition of the problem
    properties  
        % Numerical solver
        PolOrder;       % Polynomial order for each state dimension
        NumNodes;       % Number of nodes in the independent variable grid
        Basis;          % Polynomial basis 
        Grid;           % Define the independent variable grid to be used

        % State dynamics
        initial;        % Initial boundary conditions 
        final;          % Final boundary conditions
        DerDeg;         % Degree of the derivative              
        StateDim;       % State dimension 
        ControlDim;     % Control dimension

        % General parameters 
        Params;         % General parameters

        % Functions 
        BoundaryConditions;
        ControlFunction;
        CostFunction;
        LinConstraints;
        NlinConstraints; 
        InitialGuess;
        Bounds;
    end

    methods 
        % Constructor 
        function [obj] = Problem()
        end

        % Define the numerical solver to be used in the optimization process
        function [obj] = DefineSolver(obj, myOrder, myBasis, myNumNodes, myGrid)
            obj.PolOrder = myOrder; 
            obj.Basis = myBasis;
            obj.NumNodes = myNumNodes; 
            obj.Grid = myGrid; 
        end

        % Add dynamics 
        function [obj] = AddDynamics(obj, myStateDimension, myControlDimension, myL)
            obj.DerDeg = myL;
            obj.StateDim = myStateDimension;         % State dimension 
            obj.ControlDim = myControlDimension;     % Control dimension
        end

        % Add boundary conditions
        function [obj] = AddBoundaryConditions(obj, myInitial, myFinal)
            obj.initial = myInitial; 
            obj.final = myFinal;
        end

        % Add parameters 
        function [obj] = AddParameters(obj, myParameteres)
            obj.Params = myParameteres;
        end

        % Add functions 
        function [obj] = AddFunctions(obj, myBoundaryConditions, myControlFunction, myCostFunction, myLinConstraints, myNlinConstraints, myBoundsFunction, myInitialGuess)
            obj.BoundaryConditions = myBoundaryConditions;
            obj.ControlFunction = myControlFunction;
            obj.CostFunction = myCostFunction;
            obj.LinConstraints = myLinConstraints;
            obj.NlinConstraints = myNlinConstraints;
            obj.InitialGuess = myInitialGuess;
            obj.Bounds = myBoundsFunction;
        end

        % Check function 
        function [obj] = Check(obj)
            % Check the numerical solver 
            if (obj.NumNodes < 2 * max(obj.PolOrder)+1)
                warning('Quadrature may not be accurate. Consider increasing the number of nodes in the grid.');
            end

            switch (obj.Basis)
                case 'Bernstein'
                    switch (obj.Grid)
                        case 'Chebyshev'
                            error('Selected grid distribution is compatible with Bernstein polynomials');
                        case 'Legendre'
                            error('Selected grid distribution is compatible with Bernstein polynomials');
                        otherwise
                    end
 
                case 'Orthogonal Bernstein'
                    switch (obj.Grid)
                        case 'Chebyshev'
                            error('Selected grid distribution is compatible with orthogonal Bernstein polynomials.');
                        case 'Legendre'
                            error('Selected grid distribution is compatible with orthogonal Bernstein polynomials.');
                        otherwise
                    end

                case 'Chebyshev'
                    switch (obj.Grid)
                        case 'Chebyshev'
                        otherwise
                            warning('Consider selecting a Clenshaw-Curtis independent grid for maximum accuracy.');
                    end

                case 'Legendre'
                    switch (obj.Grid)
                        case 'Legendre'
                        otherwise
                            warning('Consider selecting a Lagrange-Gauss-Lobatto independent grid for maximum accuracy.');
                    end

                otherwise
            end

            % Check the dimensionality of the dynamics 
            if (size(obj.initial,1) ~= obj.StateDim * (obj.DerDeg) || size(obj.final,1) ~= obj.StateDim * (obj.DerDeg)) 
                error('Supplied boundary conditions are not of Cauchy type.');
            end

            if (length(obj.PolOrder) ~= obj.StateDim)
                warning('The input polynomial order vector mismatch the state dimension...'); 
                obj.PolOrder = [obj.PolOrder min(obj.PolOrder)*ones(1,obj.StateDim-length(obj.PolOrder))];
            elseif (size(obj.PolOrder,1) ~= obj.StateDim)
                obj.PolOrder = obj.PolOrder.';
            end

            % Check the order of the dynamics 
            if (obj.DerDeg < 1 || obj.DerDeg > 2)
                error('The problem dynamics are not well-modelled. ODE up to 2 order are supported.');
            end
        end
    end
end