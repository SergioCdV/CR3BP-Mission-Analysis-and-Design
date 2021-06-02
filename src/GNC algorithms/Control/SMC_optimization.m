%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 22/05/21
% File: SMC_optimization.m 
% Issue: 0 
% Validated: 22/05/21

%% Optimal Sliding Mode Control %%
% This script contains the function to compute the optimal parameters of an SMC controller using 
% genetic algorithms.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - string cost_function, the type of Lp norm of the control law to
%           minimize
%         - vector s0, the initial conditions of both the target and the
%           relative spacecraft
%         - scalar TOF, the time of flight

% Output: - vector parameters, the SMC controller tuned parameters

% New versions: 

function [parameters] = SMC_optimization(mu, cost_function, s0, TOF)
    %General setup 
    dt = 1e-3;              %Time step 
    tspan = 0:dt:TOF;       %Integration time span
    
    dof = 3;                %Number of parameters to tune
    PopSize = 20;          %Number of individuals at each population 
    MaxGenerations = 5;     %Maximum number of generations 
    
    %Linear constraints 
    A = []; 
    b = []; 
    Aeq = []; 
    beq = [];
    
    %Upper and lower bounds
    lb = [1e-2 0 1e-2];       %Lower bound
    ub = [1 1 0.1];           %Upper bound
    
    %General set up
    options = optimoptions(@ga,'PopulationSize', PopSize, 'MaxGenerations', MaxGenerations);

    %Compute the commands
    parameters = ga(@(param)costfunc(mu, cost_function, param, tspan, s0), dof, A, b, Aeq, beq, lb, ub, nonlcon, options);
end

%% Auxiliary functions 
%Cost function to minimize 
function [cost] = costfunc(mu, cost_function, parameters, tspan, s0)
    %Integrate the trajectory 
    GNC.Algorithms.Guidance = '';                   %Guidance algorithm
    GNC.Algorithms.Navigation = '';                 %Navigation algorithm
    GNC.Algorithms.Control = 'SMC';                 %Control algorithm
    GNC.Guidance.Dimension = 9;                     %Dimension of the guidance law
    GNC.Control.Dimension = 3;                      %Dimension of the control law
    GNC.System.mu = mu;                             %System reduced gravitational paramete
    GNC.Control.SMC.Parameters = [1 parameters];    %Controller parameters
    
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  
    [~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
    
    %Obtain the control law 
    [~, ~, u] = GNC_handler(GNC, S(:,1:6), S(:,7:12));  
    
    %Compute the control effort metrics 
    energy = zeros(3,2);                                       %Energy vector preallocation
    for i = 1:size(u,1)
        energy(i,1) = trapz(tspan, u(i,:).^2);                 %L2 integral of the control
        energy(i,2) = trapz(tspan, sum(abs(u(i,:)),1));        %L1 integral of the control
    end
    
    %Add the rendezvous error 
    lambda = 10;                                               %Lagrange multiplier
    cost = lambda*norm(S(end,7:12));
    
    switch (cost_function)
        case 'L2'
            cost = cost + norm(energy(:,1));                   %Minimize the L2 norm of the control law
        case 'L1'
            cost = cost + norm(energy(:,2));                   %Minimize the L1 norm of the control law
        otherwise
            error('No valid cost function was selected');
    end
end

%Nonlinear constraints function 
function [c, ceq] = nonlcon()
    c = []; 
    ceq = [];
end