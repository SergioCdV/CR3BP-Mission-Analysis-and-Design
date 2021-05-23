%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/05/21
% File: TOF_guidance.m 
% Issue: 0 
% Validated: 24/05/21

%% Time of flight Guidance %%
% This script contains the function to compute the optimal time of flight for a certain rendezvous maneuver.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar TOF, the time of flight for the rendezvous condition
%         - vector s0, initial conditions of both the target and the
%           relative particle
%         - scalar tol, the differential corrector scheme absolute
%           tolerance
%         - string cost_function, for both position, velocity and complete
%           rendezvous: 'Position', 'Velocity', 'State'
%         - vector sd, the desired state of the system at end of the
%           maneuver
%         - boolean two_impulsive, defining the two-impulsive/one-impulsive
%           strategy
%         - structure penalties, with the controller penalty matrices
%         - structure target_points, defining the target points for station
%           keeping
%         - structure thruster_model, defining the thruster's errro model

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [TOF] = TOF_guidance(GNC, method)


end

%% Auxiliary functions
%Nonlinear optimization core function
function [TOF] = nopt_core(cost_function, objective_norm, GNC, method)
    %Linear constraints 
    A = []; 
    b = []; 
    Aeq = []; 
    beq = [];
    
    %Upper and lower bounds
    lb = 0;                        %Lower bound
    ub = maxTOF;                   %Upper bound
    
    switch (method)
        case 'Genetic algorithm'
            %General set up
            dof = length(ub);      %Optimize the TOF
            PopSize = 100;         %Population size for each generation
            MaxGenerations = 10;   %Maximum number of generations for the evolutionary algorithm
            
            options = optimoptions(@ga,'PopulationSize', PopSize, 'MaxGenerations', MaxGenerations);
                            
            %Compute the commands
            TOF = ga(@(TOF)costfunc(cost_function, objective_norm, TOF, GNC), dof, A, b, Aeq, beq, lb, ub, ...
                     @(TOF)nonlcon(TOF,GNC), options);
            
        case 'NPL'
            %Initial guess
            sol0 = maxTOF;        
            
            %Compute the commands
            TOF = fmincon(@(TOF)costfunc(cost_function, objective_norm, TOF, GNC), sol0, A, b, Aeq, beq, lb, ub, ...
                          @(TOF)nonlcon(TOF,GNC));
            
        otherwise 
            error('No valid method was chosen');
    end
end

%Cost function 
function [cost] = costfunc(cost_function, objective_norm, TOF, GNC)
    %Switch the cost function to minimize 
    switch (cost_function)
        case 'Time'
            cost = TOF;
        case 'Control effort'
            %Compute the control law
            [u] = 1;
            
            %Control integrals
            energy = zeros(3,2);                                      %Energy vector preallocation
            for i = 1:size(u,1)
                energy(i,1) = trapz(tspan, u(i,:).^2);                %L2 integral of the control
                energy(i,2) = trapz(tspan, sum(abs(u(i,:)),1));       %L1 integral of the control
            end
            
            %Switch the Lp objective norm 
            switch (objective_norm)
                case 'L1'
                    cost = norm(energy(:,2));                         %L1 penalty
                case 'L2'
                    cost = norm(energy(:,1));                         %L2 penalty
                otherwise 
                    error('No valid Lp norm was selected');
            end
        otherwise 
            error('No valid cost function was selected');
    end
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