%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/05/21
% File: apf_guidance.m 
% Issue: 0 
% Validated: 08/05/21

%% Artificial Potential Functions Control %%
% This script contains the function to compute the control law by means of APFs.

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - structure safe_corridor, to activate/desactivate the safety APF
%         - structure Penalties, the controller scheme penalties
%         - array sO, indicating the position of the obstacles (if any) to avoid
%         - scalar TOF, the time of rendezvous to be simulated
%         - vector s0, the initial conditions of the desired guidance law

% Output: - array Sc, the rendezvous relative trajectory
%         - array dV, containing  column-wise the required number of impulses
%         - structure state, with some metrics about the scheme performance

% New versions: 

function [Sg, dV, state] = APF_control(mu, safe_corridor, Penalties, sO, TOF, s0)
    %Sanity check on initial conditions 
    if (size(s0,1) ~= 1)
        s0 = s0.';
    end
    
    %Constants 
    m = 6;                  %Phase space dimension 
    tol = 0;                %Derivative of the APF tolerance
    k = Penalties.Gain;     %Maneuver gain
    
    %Integration options 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  
    
    %Integration time 
    dt = 1e-3;                                   %Time step 
    tspan = 0:dt:TOF;                            %Integration time span
    time_horizon = length(tspan)-1;              %Final iteration time 
    
    %Preallocation 
    dV = zeros(3,length(tspan));                  %Control law
    Sg = zeros(length(tspan), 2*m);              %Controlled phase space trajectory
    
    %Initial conditions 
    Sg(1,:) = s0(1,:);

    %Integrate the desired trajectory and the control law
    for i = 1:time_horizon
        %Compute the APF evolution
        V = APF_field(safe_corridor, Penalties, sO, Sg(i,7:9).');  
        dPhi = V.'*Sg(i,10:12).';
        
        %Compute the control law
        if (dPhi >= tol)
            dV(:,i) = -k*V-Sg(i,10:12).';
        end
        
        %Change in velocity 
        Sg(i,10:12) = Sg(i,10:12) + dV(:,i).'; 
       
        %Integrate the new integration step
        [~, Saux] = ode45(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), [0 dt], Sg(i,:), options);   
        
        %Update the trajectory 
        Sg(i+1,:) = Saux(end,:);
    end
    
    %Final impulse 
    dV(:,end) = -Sg(end,10:12).';
    Sg(end,10:12) = zeros(1,3);
    
    %State output 
    state.State = true;
    state.Error = norm(Sg(end,7:12));
end

%% Auxiliary function 
% APF gradient function
function [dPhi] = APF_field(safe_corridor, Penalties, sO, S)    
    %Guidance parameters
    Q = Penalties.AttractivePenalty;             %Penalty on the distance to the origin
    R = Penalties.RepulsivePenalty;              %Penalty on the distance to the obstacles 
    sigma = Penalties.RepulsiveWidth;            %Width of the repulsive function
    
    %Compute the gradient of the attractive APF
    dPhi = Q*S;                                  %Attractive steady APF
    
    %Compute the gradient of the repulsive APF
    for i = 1:size(sO,2)
          dPhi = dPhi + exp(-((S-sO(:,i)).'*R*(S-sO(:,i)))/sigma)*(2*S-(1/sigma)*(S.'*S)*R*(S-sO(:,i)));
    end
    
    %Compute the gradient of the safety APF 
    if (safe_corridor.Safety)   
        %Safety corridor parameters
        chi = safe_corridor.Parameters(1);             %Safety corridor angle
        rho = safe_corridor.Parameters(2);             %Safety distance to the docking port
        K = safe_corridor.Parameters(3:4);             %Dimensions of the safety corridor

        %APF function
        f = S(2)^2+S(3)^2-((S(1)^3*tan(chi)^2)/(2*rho-S(1)));
        phi_s = K(1)*exp(f/K(2)); 
        dPhi = dPhi + (2*K(1)/K(2))*phi_s*[-(S(1)^2*(3*rho-S(1))*tan(chi)^2)/(2*rho-S(1))^2; S(2); S(3)];          
    else
        dPhi = dPhi + zeros(3,1); 
    end
end