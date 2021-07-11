%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/05/21
% File: MPC_control.m 
% Issue: 0 
% Validated: 20/05/21

%% Model Predictive Control Scheme %%
% This script contains the function to compute the optimal feedback control law by means of the OPTI guidance core
% and a MPC scheme.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - string cost_function, for both position, velocity and complete
%           rendezvous: 'Position', 'Velocity', 'State
%         - scalar Tmin, minimum available thrust
%         - scalar Tmax, maximum available thrust
%         - scalar TOF, the time of flight for the rendezvous condition
%         - vector s0, initial conditions of both the target and the
%           relative particle
%         - string core, selecting the solver (linear or nonlinear) to be
%           used
%         - string method, selecting the nonlinear solver to use

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [Sg, dV, state] = MPC_control(mu, cost_function, Tmin, Tmax, TOF, s0, core, method)
    %Constants 
    m = 6;                                  %Phase space dimension
    
    %Sanity check on the dimension 
    if (size(s0,1) ~= 2*m)
        s0 = s0.';                          %Initial conditions
    end
        
    %Integration setup 
    RelTol = 2.25e-14; 
    AbsTol = 1e-22; 
    options = odeset('RelTol', RelTol, 'AbsTol', AbsTol);
    
    %Initial integration    
    Phi = eye(m);                           %Initial STM
    Phi = reshape(Phi, [m^2 1]);            %Reshape the initial STM
    s0 = [s0; Phi];                         %Complete initial conditions
    dt = 1e-3;                              %Time step
    tspan = 0:dt:TOF;                       %Integration time span
    time_horizon = length(tspan)-1;         %Full MPC scheme
    
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
    St = Sn;
    Sg = St(1,:);
    
    %Preallocation 
    dV = zeros(3,length(tspan));            %Final impulsive commands
    
    %Main computation 
    for i = 1:time_horizon
        %Shrink the time span 
        atime = tspan(i:end); 
        
        %Compute the commands
        commands = OPTI_guidance(cost_function, Tmin, Tmax, atime, St, core, method);
                 
        %Add the maneuver
        dV(:,i) = shiftdim(commands(:,1));
        Sg(i,10:12) = Sg(i,10:12) + dV(:,i).';  

        %New integration
        [~, St] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), atime, Sg(i,:), options);

        %Next initial conditions
        Sg(i+1,:) = St(2,:);
    end
    
    %Final impulse
    dV(:,end) = -Sg(end,10:12);
    Sg(end,10:12) = zeros(1,3);
       
    %Final error 
    switch (cost_function)
        case 'Position'
            e = norm(St(end,7:9)); 
        case 'Velocity'
            e = norm(St(end,10:12)); 
        case 'State'
            e = norm(St(end,7:12));
        otherwise
            error('No valid cost function was selected');
    end

    state.Error = e;          %Final error 
end