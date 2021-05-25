%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 10/05/21
% File: MISS_control.m 
% Issue: 0 
% Validated: 10/05/21

%% Two-impulsive, multiple shooting control %%
% This script contains the function to compute the control law by means of an TIMS controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar TOF, the time of flight for the rendezvous condition
%         - array seed, initial orbit
%         - scalar tol, the differential corrector scheme absolute
%           tolerance
%         - number of nodes to compute the trajectory
%         - string cost_function, for both position, velocity and complete
%           rendezvous: 'Position', 'Velocity', 'State
%         - boolean two_impulsive, true if two impulses are allowed (one for
%           targetting, the second for docking)

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [Sc, dV, state] = TIMS_control(mu, TOF, seed, tol, nodes, cost_function, two_impulsive)
    %Constants 
    m = 12;                         %Complete phase space dimension 
    n = 6;                          %Individual phase space dimension
    
    %Sanity check on initial conditions dimension 
    %Sanity check on the initial seed dimensions
    if (size(seed,2) == m) || (size(seed,1) == m)
        if (size(seed,2) == m)
            seed = seed.';          %Accomodate new format
        end
    else
        error('No valid initial conditions were input');
    end
    
    %Sanity check on the cost function and the allowed number of impulses
    if (~two_impulsive)
        cost_function = 'State';
    end
    
    Phi = eye(n);                           %Initial STM  
    Phi = reshape(Phi, [n^2 1]);            %Initial STM 
    dt = 1e-3;                              %Integration time step
    h = fix(size(seed,2)/nodes)-1;          %Temporal index step
    Dt = TOF/nodes;                         %Time step between arcs
    
    switch (cost_function)
        case 'Position'
            constraints = 3;                %Additional constraints to continuity (targeting rendezvous)
        case 'Velocity'
            constraints = 3;                %Additional constraints to continuity (relative position rendezvous)
        otherwise
            constraints = 6;                %Additional constraints to continuity (full rendezvous)
    end
    
    %Preallocate internal patch points seeds 
    internalSeed = zeros((m+1)*nodes-1,1);        
    
    %Divide the orbit into the internal nodes
    for i = 1:nodes
        internalSeed(m*(i-1)+1:m*i) = seed(1:m,(i-1)*h+1);
        if (i ~= nodes)
            internalSeed(end-(nodes-1)+i) = Dt;
        end
    end   
    
    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    direction = 1;                                              %Forward integration
        
    %Differential corrector set up
    maxIter = 100;                       %Maximum number of iterations
    GoOn = true;                         %Convergence boolean 
    iter = 1;                            %Initial iteration 
    
    %Preallocation 
    ds0 = zeros(size(internalSeed,1),maxIter);                  %Vector containing the initial conditions correction
    e = zeros(m*(nodes-1)+constraints,1);                       %Error vector  
    A = zeros(m*(nodes-1)+constraints, m*nodes);                %STM matrix
    B = zeros(m*(nodes-1)+constraints, nodes-1);                %Dynamics matrix
    
    %First impulse: targeting
    while (GoOn) && (iter < maxIter)        
        for i = 1:nodes
            %Proceed with the integration
            if (i ~= nodes)
                tspan = 0:dt:internalSeed(end-(nodes-1)+i); 
            else
                tspan = 0:dt:Dt;
            end          
            S0 = shiftdim(internalSeed(m*(i-1)+1:m*i));
            S0 = [S0(1:n); Phi; S0(n+1:end); Phi];
            [~, S] = ode113(@(t,s)nlr_model(mu, direction, true, true, 'Encke', t, s), tspan, S0, options);   %New trajectory
            F = nlr_model(mu, direction, false, false, 'Encke', 0, [S(end,1:n) S(end,n+n^2+1:2*n+n^2)].');    %Vector field
            %Build the covariance matrix                                       
            if (i ~= nodes)
                %Continuity constraint
                STM = [reshape(S(end,n+1:n+n^2),[n n]) zeros(n,n); ...
                       zeros(n,n) reshape(S(end,2*n+n^2+1:end),[n n])];                 %Subarc STM
                A(m*(i-1)+1:m*i,m*(i-1)+1:m*i) = STM;                                   %Subarc STM
                A(m*(i-1)+1:m*i,m*i+1:m*(i+1)) = -eye(m);                               %Continuity constraint matrix
                B(m*(i-1)+1:m*i,i) = F(1:end);                                          %Dynamics matrix
                
                %Compute the continuity error
                e(m*(i-1)+1:m*i) = shiftdim([S(end,1:n) S(end,n+n^2+1:2*n+n^2)].'-internalSeed(m*i+1:m*(i+1)));  
            else
                %Rendezvous error
                STM = reshape(S(end,2*n+n^2+1:end),[n n]);                               %Corresponding STM
                switch (cost_function)
                    case 'Position'
                        e(end-constraints+1:end) = shiftdim(S(end,n+n^2+1:n+n^2+3)).';   %Relative phase space vector to the origin
                        A(end-constraints+1:end,end-constraints+1:end) = STM(1:3,4:6);   %Constraint matrix
                    case 'Velocity'
                        e(end-constraints+1:end) = shiftdim(S(end,n+n^2+4:2*n+n^2)).';   %Relative phase space vector to the origin
                        A(end-constraints+1:end,end-constraints+1:end) = STM(4:6,4:6);   %Constraint matrix
                    case 'State'
                        e(end-n+1:end) = shiftdim(S(end,n+n^2+1:2*n+n^2)).';             %Relative phase space vector to the origin
                        A(end-constraints+1:end,end-2:end) = STM(:,4:6);                 %Constraint matrix
                    otherwise 
                        error('No valid cost function was chosen');
                end
            end     
        end
        
        %Full covariance matrix 
        C = [A B];
                
        %Compute the correction 
        ds0(:,iter) = -pinv(C)*e;                        %Compute the variation (under-determined case)
        
        %Convergence analysis 
        if (norm(e) < tol)
            GoOn = false;
        else
            internalSeed = internalSeed+ds0(:,iter);    %Update initial conditions
            iter = iter+1;                              %Update iteration
        end       
    end
        
    %Integrate the whole trayectory
    tspan = 0:dt:sum(internalSeed(end-nodes+1:end))+Dt; 
    s0 = shiftdim(internalSeed(1:m));
    [~, S] = ode113(@(t,s)nlr_model(mu, direction, false, false, 'Encke', t, s), tspan, s0, options);
    
    %Final initial impulse 
    dV = zeros(3,length(tspan));                        %Impulses array
    dV(:,1) = seed(10:12,1)-internalSeed(n+4:m);        %Sum up each iteration contribution
    
    %Final impulse
    if (two_impulsive)
        dV(:,end) = -S(end,10:12).';    %Final impulse
        S(end,10:12) = zeros(1,3);      %Final conditions
    end

    %Output
    Sc = S;                             %Control trajectory 
    state.State = ~GoOn;                %Convergence boolean
    state.Iterations = iter;            %Number of required iterations
    state.Error = norm(e);              %Final error
end