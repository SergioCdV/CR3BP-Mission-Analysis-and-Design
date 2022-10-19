%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/05/21
% File: OPTI_guidance.m 
% Issue: 0 
% Validated: 24/05/21

%% Optimal Impulsive Guidance %%
% This script contains the function to compute the optimal control law by means of the OPTI guidance core.

% Inputs: - string cost_function, for both position, velocity and complete
%           rendezvous: 'Position', 'Velocity', 'State
%         - scalar Tmin, minimum available thrust
%         - scalar Tmax, maximum available thrust
%         - vector tspan, the time integration span
%         - array St, the trajectory along which the control law must be
%           computed
%         - string core, selecting the solver (linear or nonlinear) to be
%           used
%         - string method, selecting the nonlinear solver to use

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [commands] = OPTI_guidance(cost_function, Tmin, Tmax, tspan, St, core, method)
    %Compute the commands
    switch (core)
        case 'Nonlinear'
            commands = nopt_core(cost_function, Tmin, Tmax, tspan, St, method); 
        case 'Linear'
            commands = lopt_core(cost_function, Tmin, Tmax, St); 
        case 'Corrector'
            commands = corrector_core(cost_function, Tmin, Tmax, tspan, St);
        otherwise
            error('No valid solver was chosen')
    end
end

%% Auxiliary functions
%Differential corrector core 
function [commands] = corrector_core(cost_function, Tmin, Tmax, tspan, trajectory)
    % Differential corrector setup 
    GoOn = true; 
    iter = 1; 
    maxIter = 100; 
    tol = 1e-5; 

    mu = 0.0121505;

    % Integration tolerances (ode113)
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 

    dV = zeros(3,size(trajectory,1));

    while (GoOn && iter < maxIter)
        % Backward pass 
        [commands] = lopt_core(cost_function, Tmin, Tmax, trajectory); 

        % Forward pass 
        for i = 1:length(tspan)-1
            s0 = trajectory(i,:); 
            s0(10:12) = s0(10:12)+commands(:,i).';
           [~, aux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), [0 tspan(2)-tspan(1)], s0, options);
           trajectory(i+1,:) = aux(end,:);
        end

        e = max(sqrt(dot(dV-commands,dV-commands,1)))
        if (e < tol)
            GoOn = false; 
        else
            dV = commands;
            iter = iter+1; 
        end
    end
end

%Nonlinear optimization core function
function [commands] = nopt_core(cost_function, Tmin, Tmax, tspan, trajectory, method)
    %Constants
    impulses = length(tspan);             %Number of commands 
    
    %Linear constraints 
    A = []; 
    b = []; 
    Aeq = []; 
    beq = [];
    
    %Upper and lower bounds
    lb = Tmin*ones(1,3*impulses);         %Lower bound
    ub = Tmax*ones(1,3*impulses);         %Upper bound
    
    switch (method)
        case 'Genetic algorithm'
            %General set up
            dof = length(ub);      %Three-dimensional control vector for each instant
            PopSize = 100;         %Population size for each generation
            MaxGenerations = 10;   %Maximum number of generations for the evolutionary algorithm
            
            options = optimoptions(@ga,'PopulationSize', PopSize, 'MaxGenerations', MaxGenerations);
                            
            %Compute the commands
            solution = ga(@(u)costfunc(impulses,u), dof, A, b, Aeq, beq, lb, ub, ...
                          @(u)nonlcon(cost_function, impulses, trajectory, u), options);
                      
            commands = reshape(solution, [3 impulses]);      %Control law
            
        case 'NLP'
            %Initial guess
            sol0 = Tmax*ones(1,3*impulses);        
            
            %Compute the commands
            solution = fmincon(@(u)costfunc(impulses,u), sol0, A, b, Aeq, beq, lb, ub, ...
                               @(u)nonlcon(cost_function, impulses, trajectory, u));
            
            commands = reshape(solution, [3 impulses]);      %Control law
            
        otherwise 
            error('No valid method was chosen');
    end
end

%Linear optimization core function
function [commands] = lopt_core(cost_function, Tmin, Tmax, trajectory) 
    %Constants 
    m = 6;                                                            %Phase space dimension
    dim = size(trajectory,1);                                         %Number of impulses to make 
    Monodromy = reshape(trajectory(dim,2*m+1:end), [m m]);            %Final STM
    
    %Linear equality constraints 
    beq = -trajectory(end,m+1:2*m).';                                 %Error to rendezvous
    
    switch (cost_function)
        case 'Position'
            Aeq = zeros(m/2, m*dim);                                  %Preallocate the equality linear constraint matrix
            for i = 1:dim
                STM = reshape(trajectory(i,2*m+1:end), [m m]);        %STM at each point
                STM = Monodromy*STM^(-1);                             %Relative STM
                Aeq(:,1+(m/2)*(i-1):(m/2)*i) = STM(1:3,4:6);          %Final equality linear constraint matrix   
            end
            beq = beq(1:3);                                           %Target position
            
        case 'Velocity'
            Aeq = zeros(m/2, m*dim);                                  %Preallocate the equality linear constraint matrix
            for i = 1:dim
                STM = reshape(trajectory(i,2*m+1:end), [m m]);        %STM at each point
                STM = Monodromy*STM^(-1);                             %Relative STM
                Aeq(:,1+(m/2)*(i-1):(m/2)*i) = STM(4:6,4:6);          %Final equality linear constraint matrix   
            end
            beq = beq(4:6);                                           %Target velocity
            
        case 'State'
            Aeq = zeros(m, m*dim);                                    %Preallocate the equality linear constraint matrix
            for i = 1:dim
                STM = reshape(trajectory(i,2*m+1:end), [m m]);        %STM at each point
                STM = Monodromy*STM^(-1);                             %Relative STM
                Aeq(:,1+(m/2)*(i-1):(m/2)*i) = STM(:,4:6);            %Final equality linear constraint matrix   
            end
            
        otherwise
            error('No valid cost function was chosen');
    end
        
    %Linear inequality constraints
    A = [eye(3*dim) zeros(3*dim);  zeros(3*dim), -eye(3*dim)];
    A(1:3*dim,3*dim+1:end) = -eye(3*dim);
    A(3*dim+1:end,1:3*dim) = -eye(3*dim);
    b = zeros(m*dim,1);
    
    %Upper and lower bounds
    lb = [Tmin*ones(3*dim,1); zeros(3*dim,1)];             %Lower bound
    ub = [Tmax*ones(3*dim,1); abs(Tmax)*ones(3*dim,1)];    %Upper bound
    
    %Cost function 
    f = [zeros(1,3*dim) ones(1,3*dim)];
    
    %Solve the problem 
    options = optimset('Display', 'off');
    sol = linprog(f,A,b,Aeq,beq,lb,ub,options);
    
    %Output 
    if (isempty(sol))
        commands = zeros(3, dim);
    else
        commands = reshape(sol(1:3*dim,1), [3 dim]);
    end
end

%Cost function 
function [cost] = costfunc(impulses,u)
    %Reshape the control law 
    u = reshape(u, [3 impulses]);
    
    %Cost function 
    cost = sum(sqrt(dot(u,u,1)));
end

%Nonlinear constraints
function [c, ceq] = nonlcon(cost_function, impulses, trajectory, x)
    %Constants 
    m = 6;                  %Phase space dimension
    
    %Natural STM map
    Monodromy = reshape(trajectory(end,2*m+1:end), [m m]);          %Final STM
    e = -trajectory(end,m+1:2*m).';                                 %Final state, also the error to rendezvous
    
    %Generate the inequality function
    control = zeros(m,1);                                           %Total control effort
    for i = 1:impulses
        u = reshape(x(1,1+3*(i-1):3*i), [3 1]);                     %Control law
        STM = reshape(trajectory(i,2*m+1:end), [m m]);              %STM from the initial time to the time ti
        STM = Monodromy*STM^(-1);                                   %Relative STM from time ti to TOF
        control = control + STM(:,4:6)*u;                           %Accumulated control effort                       
    end
    
    %Nonlinear constraints
    ceq = control+e;                        %Equality constraint
    c = [];                                 %Empty inequality constraint
    switch (cost_function)
        case 'Position'
            ceq = ceq(1:3);                 %Target a position rendezvous
        case 'Velocity'
            ceq = ceq(4:6);                 %Target a velocity rendezvous
        case 'State'

        otherwise 
            error('No valid cost function was chosen');
    end
end