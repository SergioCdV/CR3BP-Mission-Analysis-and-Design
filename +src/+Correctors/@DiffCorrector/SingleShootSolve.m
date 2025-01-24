%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/01/25
% File: differential_correction.m 
% Issue: 0 
% Validated: 

%% Differential correction %%
% This function contains the implementation of the single correction
% differential corrector

% Inputs: - InitialGuess, the original orbit object to be corrcted

% Outputs: - FinalOrbit, the converged orbit object
%          - Stats, the structure containing the results of the method
 
function [FinalOrbit, Stats] = SingleShootSolve(obj, InitialGuess)
    % Constants 
    m = obj.DoFDim;             % Number of constraints in the problem

    % Newton-Rhapson method
    GoOn = true;                % Convergence flag
    iter = 1;                   % Initial iteration
    prev_e = Inf;               % Initialization of the error
    ds = zeros(m, 1);           % Initial state correction

    while (GoOn && iter < obj.Config.MaxIter)
        % Compute the error and the linear matrix
        [e, M, FinalOrbit] = obj.ComputeConFunction(InitialGuess, ds);

        % Compute the correction to the free variable
        ds = -1.0 * pinv(M) * e;

        % Convergence analysis 
        rel_convergence = norm(e - prev_e) / norm(e) <= obj.Config.RelTol;
        abs_convergence = norm(e) <= obj.Config.AbsTol;

        if ( abs_convergence || rel_convergence )
            GoOn = false;                               % Convergence is achieved
        else
            prev_e = e;                                 % Update the error
            InitialGuess = FinalOrbit;                  % Updated initial conditions
            iter = iter + 1;                            % Update iteration
        end 
    end

    if ( ~GoOn )
        % Augment initial conditions with the initial STM 
        STM = src.STM( FinalOrbit.StateDim );
        STM.Phi = eye( FinalOrbit.StateDim );
        
        % Create the complete variational system 
        VarSystem = src.Systems.VariationalCR3BP( FinalOrbit.StateDim );
        CompleteSystem = FinalOrbit.System .* VarSystem;
        
        % Integrator 
        options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14 );
        integrator = src.DynamicalSystems.HybridSolver( @ode113, options );
        
        % Initial Value Problem 
        s0 = [FinalOrbit.State(:,1); reshape(STM.Phi, [], 1)];
        VarCR3BPIVP = src.DynamicalSystems.IVP( CompleteSystem, s0, FinalOrbit.t );
        
        % Configuration 
        Solver = integrator.configure( VarCR3BPIVP );
        
        % Solve the system 
        tspan = [FinalOrbit.t(1) FinalOrbit.Period 0.01];
        [t, ~, y, ~] = Solver.solve( tspan );

        trajectory{1} = t; 
        trajectory{2} = y(1:FinalOrbit.StateDim,:);
        
        STM.Phi = y(FinalOrbit.StateDim+1:end,end);          % Monodromy matrix
        FinalOrbit.STM = STM;                                % STM of the system
        FinalOrbit.State = trajectory;                       % Final trajectory
    end

    % Final iterations 
    Stats.Iterations = iter; 
    Stats.Convergence = ~GoOn;
    Stats.RelError = norm(e - prev_e) / norm(e);
    Stats.AbsError = norm(e);
end