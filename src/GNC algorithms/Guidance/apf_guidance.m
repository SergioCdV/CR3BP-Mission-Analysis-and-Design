%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/05/21
% File: apf_guidance.m 
% Issue: 0 
% Validated: 08/05/21

%% Artificial Potential Functions Guidance %%
% This script contains the function to compute the guidance law by means of APFs.

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - structure safe_corridor, to activate/desactivate the safety APF
%         - structure Penalties, the controller scheme penalties
%         - array sO, indicating the position of the obstacles (if any) to avoid
%         - scalar TOF, the time of rendezvous to be simulated
%         - vector s0, the initial conditions of the desired guidance law
%         - boolean offline_flag, true for the offline computation of the
%           guidance law

% Output: - array output, with the definition of the guidance law depending
%           on the offline computation flag

% New versions: 

function [Sg] = APF_guidance(safe_corridor, Penalties, sO, TOF, s0, offline_flag)
    %Constants 
    m = 6;                  %Phase space dimension
    
    if (offline_flag)
        %Integration options 
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

        %Integration time 
        dt = 1e-3;                                   %Time step 
        tspan = 0:dt:TOF;                            %Integration time span

        %Preallocation of the phase space trajectory
        Sg = zeros(length(tspan),m);

        %Integrate the desired trajectory
        [~, Sg(:,1:3)] = ode45(@(t,s)APF_dynamics(safe_corridor, Penalties, sO, s), tspan, s0(1:3), options);

        for i = 1:length(tspan)
            Sg(i,4:6) = APF_dynamics(safe_corridor, Penalties, sO, Sg(i,1:3).').';
        end

        %Regress the trajectory and the derivative of the potential function
        order = 10;                                                             %Order of the approximation
        [Cp, Cv, Cg] = CTR_guidance(order, tspan, Sg);
        Sg = [Cp; Cv; Cg]; 
    else
        Sg(1:3) = zeros(3,1);
        Sg(4:6) = APF_dynamics(safe_corridor, Penalties, sO, s0(1:3));
        Sg(7:9) = zeros(3,1);
    end
end

%% Auxiliary function 
%APF evaluation
function [Phi] = APF_evaluation(safe_corridor, Penalties, sO, S)
    %Guidance parameters
    Q = Penalties.AttractivePenalty;                   %Penalty on the distance to the origin
    R = Penalties.RepulsivePenalty;                    %Penalty on the distance to the obstacles 
    sigma = Penalties.RepulsiveWidth;                  %Width of the repulsive function
    
    %Compute the attractive APF value
    phi_a = (1/2)*S.'*Q*S;                             %Attractive steady APF

    %Compute the repulsive APF values 
    phi_r = 0;
    for i = 1:size(sO,2)
        phi_r = phi_r + (S.'*S)*exp(-((S-sO(:,i)).'*R*(S-sO(:,i)))/sigma);   
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

%APF gradient function
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

%APF dynamics 
function [dS] = APF_dynamics(safe_corridor, Penalties, sO, S)
    dS = -APF_field(safe_corridor, Penalties, sO, S(1:3));   
end