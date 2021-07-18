%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/05/21
% File: TOF_guidance.m 
% Issue: 0 
% Validated: 24/05/21

%% Time of flight Guidance %%
% This script contains the function to compute the optimal time of flight for a certain rendezvous maneuver.

% Inputs: - string cost_function, to optimize for the minimum time or the
%           time of minimum control effort
%         - string objective_norm, the Lp norm of the control law to minimize 
%         - function handle control_scheme, the controller scheme with
%           whom the TOF guidance must be evaluated
%         - scalar maxTOF, the upper limit of the rendezvous time
%         - string method, the nonlinear core to solve the optimization
%           problem

% Output: - scalar TOF, the optimal time of flight

% New versions: 

function [TOF] = TOF_guidance(cost_function, objective_norm, control_scheme, maxTOF, method)
    %Linear constraints
    A = []; 
    b = []; 
    Aeq = []; 
    beq = [];
    
    %Upper and lower bounds
    lb = 2e-3;                     %Lower bound
    ub = maxTOF;                   %Upper bound
    
    switch (method)
        case 'Genetic algorithm'
            %General set up
            dof = length(ub);      %Optimize the TOF
            PopSize = 100;         %Population size for each generation
            MaxGenerations = 10;   %Maximum number of generations for the evolutionary algorithm
            
            options = optimoptions(@ga,'PopulationSize', PopSize, 'MaxGenerations', MaxGenerations);
                            
            %Compute the commands
            TOF = ga(@(TOF)costfunc(cost_function, objective_norm, control_scheme, TOF), dof, A, b, Aeq, beq, lb, ub, ...
                     @(TOF)nonlcon(control_scheme, TOF), options);
            
        case 'NLP'
            %Initial guess
            sol0 = maxTOF;        
            
            %Compute the commands
            TOF = fmincon(@(TOF)costfunc(cost_function, objective_norm, control_scheme, TOF), sol0, A, b, Aeq, beq, lb, ub, ...
                          @(TOF)nonlcon(control_scheme, TOF));
            
        otherwise 
            error('No valid method was chosen');
    end
end

%Cost function 
function [cost] = costfunc(cost_function, objective_norm, control_scheme, TOF)
    %Regularize with the rendezvous error
    dt = 1e-3;                                  %Time step
    tspan = 0:dt:TOF;                           %Integration time
    [~, u] = feval(control_scheme, TOF);        %Compute the control law   
                    
    %Switch the cost function to minimize 
    switch (cost_function)
        case 'Time'
            cost = TOF;
            
        case 'Control effort' 
            %Control integrals
            effort = control_effort(tspan, u, true);          %Control effort figures of merit
                        
            %Switch the Lp objective norm 
            switch (objective_norm)
                case 'L1'
                    cost = norm(effort(:,2));           %L1 penalty
                case 'L2'
                    cost = norm(effort(:,1));           %L2 penalty
                otherwise 
                    error('No valid Lp norm was selected');
            end
                        
        otherwise 
            error('No valid cost function was selected');
    end
end

%Nonlinear constraints
function [c, ceq] = nonlcon(control_scheme, TOF)    
    %Constants 
    tol = 1e-5;       %Rendezvous error tolerance 
    
    %Integrate the trajectory 
    dt = 1e-3;                                  %Time step
    tspan = 0:dt:TOF;                           %Integration time
    [S, ~] = feval(control_scheme, TOF);        %Compute the control law   
    
    %Figures of merit 
    e = figures_merit(tspan, S); 
    
    %Nonlinear constraints
    ceq = [];                        %Equality constraint
    c = e(end)-tol;                  %Inequality constraint
end