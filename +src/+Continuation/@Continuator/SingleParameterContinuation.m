%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/01/25
% File: Continuator.m 
% Issue: 0 
% Validated: 

%% Single Parameter Continuator %% 
% This class implements a generic continuator object through SPC %

function [ObjectFamily] = SingleParameterContinuation(obj, InitialGuess)
    % Pre-allocation 
    ObjectFamily = {};      % Pre-allocation of the family

    % Initial step length
    ds = 0;
    [e, ~] = obj.ContinuationFunction(InitialGuess, ds);
    ds = e / obj.Config.MaxIter;

    % Main loop 
    iter = 1;               % Number of iterations of the method
    idx = 1;                % Index of the family
    GoOn = true;            % Convergence boolean
    prev_e = Inf;           % Initialization of the relative error function

    while ( GoOn && iter < obj.Config.MaxIter)
        % Update initial conditions
        [e, UpdatedObject] = obj.ContinuationFunction(InitialGuess, ds);

        % Refine the solution
        [RefinedObject, Stats] = obj.DiffCorrector( UpdatedObject );

        % Convergence analysis
        rel_convergence = norm(e - prev_e, 'inf') / norm(e, 'inf') <= obj.Config.RelTol;
        abs_convergence = norm(e) <= obj.Config.AbsTol;

        if ( abs_convergence || rel_convergence )
            GoOn = false;                               % Stop the continuation process

        elseif ( ~Stats.Convergence )
            ds = ds / 10;                               % Reduce the step of the continuation

        else
            prev_e = e;                                 % Update the error
            InitialGuess = RefinedObject;               % New initial guess 
            ds = e / (obj.Config.MaxIter - iter);       % Update the continuation step

            ObjectFamily{idx} = RefinedObject;          % Save the object
            idx = idx + 1;                              % Update the family index
        end
    end
end