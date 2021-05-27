%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/05/21
% File: apf_guidance.m 
% Issue: 0 
% Validated: 08/05/21

%% Artificial Potential Functions guidance %%
% This script contains the function to compute the guidance law by means of APFs.

% Inputs: - string dynamics, to select the dsired type of APF algorithm
%         - boolean safe_corridor, to activate/desactivate the safety APF
%         - vector state, the state of the object from which the guidance
%           law is computed
%         - array obstacles_state, indicating the position of the obstacles
%           (if any) to avoid
%         - scalar dt, the time step taken in the simulation 
%         - vector Sg0, the initial conditions of the desired guidance law

% Output: - the index s, containing information about close bifurcations around the solution associated with the STM.

% New versions: 

function [Cp, Cv, Cg, dPhi] = APF_guidance(dynamics, safe_corridor, Penalties, sO, TOF, s0)
    %Integration options 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  
    
    %Integration time 
    dt = 1e-3;                                   %Time step 
    tspan = 0:dt:TOF;                            %Integration time span
       
    %Preallocation 
    dPhi = zeros(length(tspan));          %Preallocation of the APF total derivatives

    %Integrate the desired trajectory
    [~, Sg] = ode45(@(t,s)APF_dynamics(dynamics, safe_corridor, Penalties, sO, s), tspan, s0(1:3), options);                                               

    %Compute the APF evolution
    for i = 1:length(tspan)
       dPhi(i) = APF_evaluation(dynamics, safe_corridor, Penalties, sO, Sg(i,1:3).');  
    end
            
    %Regress the trajectory and the derivative of the potential function
    order = 10;                                                             %Order of the approximation
    [Cp, Cv, Cg] = CTR_guidance(order, tspan, Sg);
end

%% Auxiliary function 
%APF evaluation
function [Phi] = APF_evaluation(dynamics, safe_corridor, Penalties, sO, S)
    %Guidance parameters
    Q = Penalties.AttractivePenalty;             %Penalty on the distance to the origin
    R = Penalties.RepulsivePenalty;              %Penalty on the distance to the obstacles 
    
    %Compute the attractive APF value
    switch (dynamics)
        case 'Steady'
            phi_a = (1/2)*S.'*Q*S;               %Attractive steady APF
        case 'Unsteady'
            error('Selected algorithm not implemented');
        otherwise
            error('No valid APF dynamics were selected');
    end

    %Compute the repulsive APF values 
    phi_r = 0;
    for i = 1:size(sO,2)
        phi_r = phi_r + (1/2)*(S.'*Q*S)/((S-sO(:,i)).'*R*(S-sO(:,i))-1);   
    end

    %Compute the safety APF value
    if (safe_corridor.Safety)   
        chi = safe_corridor.Parameters(1);             %Safety corridor angle
        rho = safe_corridor.Parameters(2);             %Safety distance to the docking port
        K = safe_corridor.Parameters(3:4);             %Dimensions of the safety corridor

        %APF function
        f = S(2)^2+S(3)^2-((S(1)^3*tan(chi)^2)/(2*rho-s(1)));
        phi_s = K(1)*exp(-abs(f)/K(2));  
    else
        phi_s = 0; 
    end
    
    %Evaluate the total APF 
    Phi = phi_a + phi_s + phi_r; 
end

%APF vector field 
function [dPhi] = APF_field(dynamics, safe_corridor, Penalties, sO, S)
    %Initialization 
    dPhi = zeros(3,1); 
    
    %Guidance parameters
    Q = Penalties.AttractivePenalty;             %Penalty on the distance to the origin
    R = Penalties.RepulsivePenalty;              %Penalty on the distance to the obstacles 
    
    %Compute the gradient of the attractive APF
    switch (dynamics)
        case 'Steady'
            dPhi = dPhi + Q*S;                 %Attractive steady APF
        case 'Unsteady'
            error('Selected algorithm not implemented');
        otherwise
            error('No valid APF dynamics were selected');
    end
    
    %Compute the gradient of the repulsive APF
    for i = 1:size(sO,2)
        dPhi = dPhi + (1/2)*(((S-sO(:,i)).'*R*(S-sO(:,i))-1)*Q*S-S.'*Q*S*R*(S-sO(:,i))) ...
                      /((S-sO(:,i)).'*R*(S-sO(:,i))-1)^2;   
    end
    
    %Compute the gradient of the safety APF 
    if (safe_corridor.Safety)   
        chi = safe_corridor.Parameters(1);             %Safety corridor angle
        rho = safe_corridor.Parameters(2);             %Safety distance to the docking port
        K = safe_corridor.Parameters(3:4);             %Dimensions of the safety corridor

        %APF function
        phi_s = K(1)*(exp(-(1/K(2))*(S(2)^2+S(3)^2-((S(1)^3*tan(chi)^2)/(S*rho-s(1))))));  
    else
        dPhi = dPhi + zeros(3,1); 
    end
end

%APF dynamics 
function [dS] = APF_dynamics(dynamics, safe_corridor, Penalties, sO, S)
    dS = -APF_field(dynamics, safe_corridor, Penalties, sO, S(1:3));   
end