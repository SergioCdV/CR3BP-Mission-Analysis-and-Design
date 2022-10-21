%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/05/21
% File: SMC_control.m 
% Issue: 0 
% Validated: 08/05/21

%% Sliding Mode Control %%
% This script contains the function to compute the control law by means of an SMC controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - array Sg, the guidance law to follow
%         - array Sn, the system state
%         - vector parameters, defining the tuned parameters of the controller

% Output: - vector u, the computed control law

% New versions: 

function [u] = SMC_control(mu, Sg, Sn, parameters, model)
    %SMC parameters 
    lambda = parameters(1);                             %General loop gain
    epsi = parameters(2);                               %Reachability condition gain
    alpha = parameters(3);                              %Reachability condition exponent
    delta = parameters(4);                              %Boundary layer width
    
    %Preallocation 
    u = zeros(3,size(Sn,1));                            %Control vector
    
    for i = 1:size(Sn,1)
        %Relative orbital state
        r = Sn(i,7:9).';                                %Instanteneous position vector
        v = Sn(i,10:12).';                              %Instanteneous velocity vector

        %Compute the position and velocity errors
        dr = r-Sg(i,1:3).';                             %Position error
        dv = v-Sg(i,4:6).';                             %Velocity error

        %Torque computation
        s = dv+lambda*dr;                               %Sliding surface

        f = nlr_model(mu, true, false, false, model, 0, Sn(i,:).');               %Relative CR3BP 
        
        %Final control law
        ds = epsi*(norm(s)^(alpha)*saturation(s, delta).'+s);                      %Reachability condition function
        u(:,i) = Sg(i,7:9).'-f(10:12)-lambda*dv-ds;                                %Control vector
    end
end

%% Auxiliary functions
%Saturation function
function [U] = saturation(s, delta)
    %Compute a bang bang saturation law, given the boundary layer delta 
    U = zeros(1,length(s));
    for i = 1:size(s)
        if (s(i) > delta)
            U(i) = 1; 
        elseif (s(i) < -delta)
            U(i) = -1; 
        else
            U(i) = (1/delta)*s(i);
        end
    end
end