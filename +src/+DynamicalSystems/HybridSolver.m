%% KSHLand %% 
% Author: Sergio Cuevas
% Date: 15/11/2024

%% Hybrid Solver %% 
% This is a class defining a generic hybrid system solver %

classdef HybridSolver
    properties
        integrator;         % Generic integrator 
        int_options;        % Integration options 
        problem;            % Problem to be solved
    end

    methods
        % Basic constructor
        function [obj] = HybridSolver( integrator_handler, int_options )
            % Sanity checks
            if ( ~isa(integrator_handler, 'function_handle') )
                error('The input integrator set cannot be evaluated... Aborting')
            else
                % Supported integrator classes
                solversOde = {'ode45', 'ode23', 'ode113', 'ode15s', 'ode23s', 'ode23t', 'ode23tb', 'ode78', 'ode89', 'ode15i'};
                
                % Check if the selected integrator belong to the supported 
                funcName = func2str(integrator_handler); 
                esOde = ismember(funcName, solversOde);

                if ( esOde )
                    obj.integrator = integrator_handler;
                else
                    error('The selected integrator is not supported... Aborting')
                end
            end
            
            obj.int_options = int_options;
        end

        % Configuration 
        function [obj] = configure(obj, myProblem )
            % Problems
            obj.problem = myProblem;
        end

        % Solver 
        function [t, x, stats] = solve(obj, tspan, jspan) 
            % Sanity checks 
            if (~exist('tspan', 'var') )
                tspan = [obj.problem.t0 1E5];
            else
                time_span = tspan(1):tspan(3):tspan(2);
                t = zeros(1, length(time_span) + 1);                       % Independent variable 
                x = zeros( obj.problem.System.StateDim, length(t) );       % State of the system
            end

            if (~exist('jspan', 'var') )
                jspan = [0 1E5];

                if (~exist('tspan', 'var') )
                    error('Either the continuous horizon or the jump horizon must be definite... Aborting');
                else
                    j = zeros( size(t) ); 
                end
            else 
                index_span = 1:jspan(end);
                j = zeros(1, length(index_span) + 1);                      % Independent counter 

                if ( isinf(tspan(end)) )
                    x = zeros( obj.problem.System.StateDim, length(j) );   % State of the system
                    t = zeros( size(j) );
                else
                    max_steps = min(length(t), length(j));
                    x = x(:,1:max_steps);                                  % State of the system
                end
            end

            % Pre-allocation 
            stats = [];

            % Initialization 
            step = 1; 
            x(:,step) = obj.problem.IC;     % Initial state conditions
            t(step) =   obj.problem.t0;     % Initial clock

            % Sanity check 
            if ( tspan(1) ~= t(1) )
                error('The initial conditions are not synchronized with the independent variable... Aborting');
            end
            
            % Solving of the problem
            while ( t(step) < tspan(2) && j(step) < jspan(end) )
                % Check for pre-allocation 
                if ( step + 1 > min( size(t,2), size(j,2) ) )
                    x = [x zeros( size(x,1), 1E5 )];
                    t = [t zeros( 1, 1E5 )];
                    j = [j zeros( 1, 1E5 )];
                end

                % Check if we are in the jump set or flow set
                controller = obj.problem.System.ExogenousInput(t(step), j(step), x(:,step), obj.problem.System.params);
                flow_flag = obj.problem.System.inFlowSet( t(step), j(step), x(:,step), controller, obj.problem.System.params );
                jump_flag = obj.problem.System.inJumpSet( t(step), j(step), x(:,step), controller, obj.problem.System.params );

                flowing = flow_flag && (obj.problem.System.PriorityRule == 2 || (obj.problem.System.PriorityRule == 1 && ~jump_flag));
                jumping = jump_flag && (obj.problem.System.PriorityRule == 1 || (obj.problem.System.PriorityRule == 2 && ~flow_flag));

                if ( flowing )
                    % Function handle for the controller 
                    controller = @(t,x)obj.problem.System.ExogenousInput(t, j(step), x, obj.problem.System.params);

                    % Prepare the integration 
                    obj.int_options = odeset(obj.int_options, 'Events', @(t, x)obj.event(t, j(step), x, controller(t, x), obj.problem.System.params));
                    
                    % Integration
                    time_step = t(step) : tspan(3) : (tspan(2) + tspan(3)); 
                    [t_aux, x_aux, ~, ~, ie] = obj.integrator( @(t,x)obj.problem.System.Dynamics(t, j(step), x, controller(t, x), obj.problem.System.params), time_step, x(:,step), obj.int_options ); 
                    x_aux = x_aux.';
                
                    % Check if the event takes at the second step 
                    controller = obj.problem.System.ExogenousInput(t_aux(2), j(step), x_aux(:,2), obj.problem.System.params);
                    missed_flow = obj.problem.System.inFlowSet( t_aux(2), j(step), x_aux(:,2), controller, obj.problem.System.params );
                    missed_jump = obj.problem.System.inJumpSet( t_aux(2), j(step), x_aux(:,2), controller, obj.problem.System.params );
                    missed_event = missed_flow && ~(missed_jump && obj.problem.System.PriorityRule == 1);

                    if ( missed_event || length(t_aux) == 2 )
                        % Save the values
                        num_steps = length(t_aux) - 1;
    
                        x(:, step + 1: step + num_steps) = x_aux(:,2:end);
                        t(step + 1: step + num_steps) = t_aux(2:end);
                        j(step + 1: step + num_steps) = j(step) * ones(1, num_steps);

                        step = step + num_steps;

                    else
                        dt = t_aux(2) - t_aux(1);
                        delta = 10^(-9:9);
                        GoOn = true; 
                        iter = 1;

                        while (GoOn)
                            if (dt <= delta(iter))
                                t_plus = t_aux(2);
                                x_plus = x_aux(2);
                                GoOn = false;
                            else
                                % Euler forward 
                                t_plus = t(step) + delta(iter);
                                controller = obj.problem.System.ExogenousInput(t_plus, j(step), x(:,step), obj.problem.System.params);
                                x_plus = x(:,step) + delta(iter) * obj.problem.System.Dynamics(t_plus, j(step), x(:,step), controller, obj.problem.System.params);
    
                                % Check if we are leaving the flow set
                                missed_flow = obj.problem.System.inFlowSet( t_plus, j(step), x_plus, controller, obj.problem.System.params );
                                missed_jump = obj.problem.System.inJumpSet( t_plus, j(step), x_plus, controller, obj.problem.System.params );
                                missed_event = missed_flow && ~(missed_jump && obj.problem.System.PriorityRule == 1);
    
                                if (missed_event)
                                    GoOn = false;
                                else
                                    iter = iter + 1;
                                end
                            end
                        end

                       % Save results
                       step = step + 1;
                       t(step) = t_plus; 
                       x(:,step) = x_plus;
                       j(step) = j(step-1);
                    end

                % Check if we are in the jump set
                elseif ( jumping )
                    % Jump
                    controller = obj.problem.System.ExogenousInput( t(step), j(step), x(:,step), obj.problem.System.params );
                    [~, xp] = obj.problem.System.Jump( t(step), j(step), x(:,step), controller, obj.problem.System.params );

                    % Save the values
                    step = step + 1; 
                    x(:,step) = xp;
                    t(step) = t(step - 1);
                    j(step) = j(step - 1) + 1;
                end
            end

            % Restrict the results 
            x = x(:, 1:step); 
            t = t(:, 1:step);
            j = j(:, 1:step);
        end
        
        % Halt events
        function [val, isterminal, direction] = event(obj, t, j, x, u, params)
            % Sets of the problem
            jump = obj.problem.System.inJumpSet(t, j, x, u, params);  
            flow = obj.problem.System.inFlowSet(t, j, x, u, params);

            switch ( obj.problem.System.PriorityRule )
                case 1
                    % For jump priority, we terminate flows whenever (1) the solution leaves C, or (2) the solution enters D
                    stop = ~flow || jump;
                case 2
                    % For flow priority, we terminate flows whenever the solution leaves C
                    stop = ~flow;
            end
            
            if any(isnan(x)) || any(isinf(x))
                stop = 1;
            end

            val = 1 - stop;     % Value used to terminate the flow
            isterminal = 1;     % Terminate the integration
            direction = -1;     % Terminate for decreasing values
        end
    end
end