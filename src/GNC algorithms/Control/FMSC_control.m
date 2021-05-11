%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 11/05/21
% File: FMSC_control.m 
% Issue: 0 
% Validated: 11/05/21

%% Floquet Mode Safe Control %%
% This script contains the function to compute the control law by means of an FMSC controller.

% Inputs: - string model, selecting the linear model to compute the linear state
%           transition matrix of the system
%         - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - array St, the target orbit reference state 
%         - array Sg, the guidance law to follow
%         - array Sn, the system state
%         - scalar Ln, the libration point number. It may be left 0 if the
%           target orbit is not librating around any Lagrange point
%         - scalar gamma, the relative distance of the Ln point to the
%           nearest primary. Again, it may be left as 0 if needed
%         - matrices Q and M, penalizing on the state error and the control
%           effort

% Output: - vector u, the computed control law

% New versions: 

function [Sc, dV, state] = FMSC_control(mu, TOC, s0, constraint, restriction)
    %Constants 
    m = 6;       %Phase space dimension
    
    %Sanity check on initial conditions dimension 
    if (size(s0,1) == 1)
        s0 = s0.';
    end
    
    %Preallocation 
    J = zeros(2,length(tspanc)-1);                              %Cost function to analyze
    dV = zeros(length(tspanc)-1,3);                             %Velocity impulses all along the look ahead time arc

    %Select the restriction level of the CAM 
    constrained = constraint.Constrained;
    lambda = constraint.SafeDistance;                           %Safety distance
    
    %Integration tolerance and setup 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22); 
    
    dt = 1e-3;                                                  %Integration time step  
    tspan = 0:dt:TOC;                                           %Integration time span
    
    %Set up the differential corrector
    maxIter = 100;                                     %Maximum number of iterations
    GoOn = true;                                       %Convergence boolean 
    iter = 1;                                          %Initial iteration 
    
    %Initial conditions and integration
    Phi = eye(m);                                      %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                       %Reshape the initial STM
    s0 = [s0; Phi];                                    %Complete phase space + linear variational initial conditions
    
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);

    for i = 1:length(tspanc)-1
        %Shrink the look ahead time 
        atime = tspan(i:end);

        %Compute the Floquet modes at each time instant 
        Monodromy = reshape(Sn(end,13:end), [m m]);                     %State transition matrix at each instant        
        [E, sigma] = eig(Monodromy);                                    %Eigenspectrum of the STM 
        Phi = Monodromy*reshape(S(i,13:end), [m m])^(-1);               %Relative STM

        for j = 1:size(E,2)
            if (constrained)
                E(:,j) = sigma(j,j)*E(:,j);
            else
                E(:,j) = exp(-tspan(i)/tspan(end)*log(sigma(j,j)))*E(:,j);
            end
        end

        %Compute the maneuver
        switch (restriction) 
            case 'Worst'
                if (constrained)
                    safeS = lambda(1)*E(:,1);                                      %Safety constraint
                    STM = Phi(:,4:6);                                              %Correction matrix
                    error = safeS-S(end,7:12).';                                   %Error in the unstable direction only
                    maneuver = pinv(STM)*error;                                    %Needed maneuver 
                    dV(:,i) = maneuver(end-2:end);                                 %Needed change in velocity
                else
                    w = [0; 1; 1; 1];                                              %Some random vector offerint triaxial control
                    STM = [E(:,1) -[zeros(3,3); eye(3)]];                          %Linear application
                    error = lambda(1)*rand(6,1);                                   %State error vector
                    maneuver = pinv(STM)*error;                                    %Needed maneuver
                    dV(:,i) = real(maneuver(end-2:end));                           %Needed change in velocity
                end
            case 'Best'
                if (constrained)
                    safeS = lambda(1)*E(:,1)+lambda(2).*E(:,3:end);                %Safety constraint
                    STM = Phi(:,4:6);                                              %Correction matrix
                    error = safeS-S(end,7:12).';                                   %Error in the unstable direction only
                    maneuver = pinv(STM)*error;                                    %Needed maneuver 
                    dV(:,i) = maneuver(end-2:end);                                 %Needed change in velocity
                else
                    STM = [E(:,1) E(:,3:end) -[zeros(3,3); eye(3)]];               %Linear application
                    error = 1e-3*rand(6,1);                                        %State error vector
                    maneuver = STM.'*(STM*STM.')^(-1)*error;                       %Needed maneuver
                    dV(:,i) = real(maneuver(end-2:end));                           %Needed change in velocity
                end
            otherwise
                error('No valid case was selected');
        end

        %Integrate the trajectory 
        s0 = Sn(i,1:12);                        %Initial conditions
        s0(10:12) = s0(10:12)+real(dV(i,:)).';  %Update initial conditions with the velocity change

        [~, s] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), atime, s0, options);
        
        %Evaluate the cost function
        J(1,i) = (1/2)*(s(end,7:9)-so(1,1:3))*Q*(s(end,7:9)-so(1,1:3)).';
        J(2,i) = norm(maneuver)+(1/2)*s(end,7:12)*s(end,7:12).';  
     end
end