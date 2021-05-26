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

function [Sg, phi] = APF_guidance(dynamics, safe_corridor, Penalties, obstacles_states, TOF, s0)
    %Constants 
    m = 6;                                       %State dimension 
            
    %Integration options 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  
    
    %Integration time 
    dt = 1e-3;                                   %Time step 
    tspan = 0:dt:TOF;                            %Integration time span
   
    %Switch between the two types of control schemes 
    switch (controller_scheme)
        case 'Impulsive'
            %Preallocation 
            Sg = zeros(length(tspan), m+3);      %Preallocation of the reference phase space trajectory
            phi = zeros(length(tspan));          %Preallocation of the APF
            
            %Integrate the desired trajectory
            [~, S] = ode45(@(t,s)APF_dyanmics(), tspan, s0, options);       %Reference position integration 
            Sg(:,1:m) = S;                                                  %Reference position and velocity
            
            %Compute the acceleration field 
            for i = 1:length(tspan)
                phi(i) = APF_evaluation();
                Sg(i,7:9) = APF_dynamics();                                  
            end
            
        case 'Continuous'        
            %Preallocation 
            Sg = zeros(length(tspan), m/2);      %Preallocation of the reference phase space trajectory
            phi = zeros(length(tspan));          %Preallocation of the APF
            
            %Compute the desired velocity
            [~, S] = ode45(@(t,s)APF_dyanmics(), tspan, s0, options);       %Reference position integration
            
            for i = 1:length(tspan)
                phi(i) = APF_evaluation();      
                Sg(i,:) = S(i,4:6);              
            end
            
        otherwise
            error('No valid controller was chosen');
    end        
end

%% Auxiliary function 
%APF evaluation
function [Phi] = APF_evaluation(dynamics, safe_corridor, Penalties, sO, S)
    %Guidance parameters
    Q = Penalties.AttractivePenalty;             %Penalty on the distance to the origin
    R = Penalties.AttractivePenalty;             %Penalty on the distance to the obstacles 
    
    %Compute the attractive APF value
    switch (dynamics)
        case 'Steady'
            phi_a = (1/2)*S*Q*S.';                              %Attractive steady APF
        case 'Unsteady'
            error('Algorithm not implemented');
        otherwise
            error('No valid APF dynamics were selected');
    end

    %Compute the repulsive APF values 
    obs = reshape(sO(1:3,:), 3*size(sO,2), 1);
    phi_r = (1/2)*(S.'*Q*S)/((repmat(S,size(sO,2),1)-obs).'*R*(repmat(S,size(sO,2),1)-obs)-1);    

    %Compute the safety APF value
    if (safe_corridor.Safety)   
        chi = safe_corridor.Parameters(1);             %Safety corridor angle
        rho = safe_corridor.Parameters(2);             %Safety distance to the docking port
        K = safe_corridor.Parameters(3:4);             %Dimensions of the safety corridor

        %APF function
        phi_s = K(1)*(exp(-(1/K(2))*(S(2)^2+S(3)^2-((S(1)^3*tan(chi)^2)/(S*rho-s(1))))));  
    else
        phi_s = 0; 
    end
    
    %Evaluate the total APF 
    Phi = phi_a + phi_s + phi_r; 
end

%APF dynamics 
function [dS] = APF_dynamics(dynamics, safe_corridor, Penalties, sO, S)    
    %Evaluate the dynamics numerically
    dS(1:3) = -S(4:6);                                          
    dS(4:6) = -numerical_jacobian(1, @(s)APF_evaluation(dynamics, safe_corridor, Penalties, sO, s), ...
                                  S(1:3));   
end