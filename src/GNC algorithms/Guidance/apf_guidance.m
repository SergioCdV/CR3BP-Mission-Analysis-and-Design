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

function [Sg, phi] = APF_guidance(dynamics, safe_corridor, Sn, obstacles_states, dt)
    %Constants 
    m = 6;                                       %State dimension 
    Q = eye(m/2);                                %Penalty on the distance to the origin
    R = eye(m/2*size(obstacles_states,2));       %Penalty on the distance to the obstacles 
        
    %Integration options 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);   
    
    %Preallocation 
    Sg = zeros(size(Sn,1), m+3);                 %Guidance law
    Sg(1,1:3) = Sn(1,1:3);                       %Initial conditions

    for i = 1:size(Sn,1)
        %Compute the attractive APF value
        switch (dynamics)
            case 'Steady'
                phi = (1/2) * Sn(i,1:3)*Q*Sn(i,1:3).';       %Attractive steady APF
                dPhi = Q*Sn(i,1:3).';                        %Gradient of the APF 
                hPhi = Q;                                    %Hessian of the APF 
            case 'Unsteady'
                error('Algorithm not implemented');
            otherwise
                error('No valid APF dynamics were selected');
        end

        %Compute the repulsive APF values 
        obs = reshape(obstacles_states(1:3,:), 3*size(obstacles_states,2), 1);
        phi_r = @(s)(1/2)*(s.'*Q*s)/((repmat(s,size(obstacles_states,2),1)-obs).'*R*(repmat(s,size(obstacles_states,2),1)-obs)-1);    

        %Compute the safety APF value
        if (safe_corridor.Safety)   
            chi = safe_corridor.Parameters(1);             %Safety corridor angle
            rho = safe_corridor.Parameters(2);             %Safety distance to the docking port
            K = safe_corridor.Parameters(3:4);             %Dimensions of the safety corridor

            %APF function
            phi_s = @(s)K(1)*(exp(-(1/K(2))*(s(2)^2+s(3)^2-((s(1)^3*tan(chi)^2)/(2*rho-s(1))))));  
        else
            phi_s = @(s)(0); 
        end

        %Compute the guidance trajectory
        Sg(i,4:6) = -(dPhi + numerical_jacobian(1, phi_r, Sn(i,1:3).').' + numerical_jacobian(1, phi_s, Sn(i,1:3).').');    %Reference velocity
        [~, r] = ode45(@(t,r)(Sg(i,4:6).'), [0 dt], Sn(i,1:3), options);                                                    %Reference position integration 
        Sg(i,1:3) = r(end,:).';                                                                                             %Reference position
        Sg(i,7:9) = (hPhi + numerical_hessian(phi_r, Sn(i,1:3).') + numerical_hessian(phi_s, Sn(i,1:3).'))*Sg(i,4:6).';     %Reference acceleration
    end
end