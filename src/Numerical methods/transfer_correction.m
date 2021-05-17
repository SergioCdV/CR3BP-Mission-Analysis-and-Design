%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 03/01/21
% File: differential_correction.m 
% Issue: 0 
% Validated: 

%% Transfer differential correction %%
% This function contains the algorithm to compute several correction
% algorithms to generate low-budget transfers in the C3RBP context.

% Inputs: - string algorithm, selecting the differential correction scheme to use. 
%         - double mu, the reduced gravitational parameter of the system.
%         - object parking_orbit, defining the initial parking orbit.
%         - object target_orbit, which vary depending on the algorithm in use. 
%         - int maxIter, number of maximum allowed corrections.
%         - double tol, tolerance to stop the correction. 

% Outputs: - vector xf, converged/last corrected trajectory as well as other related parameter such as the TOF.
%          - boolean state, outputing the convergence of the algorithm.

% Methods: 

% New versions: 

function [xf, dVf, state] = transfer_correction(algorithm, mu, parking_orbit, target_orbit, maxIter, tol, varargin)            
    %Implement the selected scheme 
    switch (algorithm)
        case 'HOI transfer'
            [xf, dVf, state] = HOI_transfer(mu, parking_orbit, target_orbit, maxIter, tol, varargin);
        otherwise
            error('No valid option was selected');
    end
end

%% Auxiliary functions 
%Compute a HOI transfer using the stable manifold, as detailed in Barden, 1994
function [xf, dVf, state] = HOI_transfer(mu, parking_orbit, target_orbit, maxIter, tol, varargin)
    %Constants 
    m = 6;                                      %Phase space dimension
    
    %Parking orbit definition 
    Primary = parking_orbit.Primary;            %Primary associated to the parking orbit   
    switch (Primary)
        case 'First'
            R = [-mu; 0; 0];                    %Primary location in the configuration space
            branch = 'L';                       %Manifold branch to globalize
            map = 'First primary';              %Poincaré map to use
            event = @(t,s)fp_crossing(t,s,mu);  %Integration event
        case 'Secondary'
            R = [1-mu; 0; 0];                   %Primary location in the configuration space
            branch = 'R';                       %Manifold branch to globalize
            map = 'Secondary primary';          %Poincaré map to use
            event = @(t,s)sp_crossing(t,s,mu);  %Integration event
        otherwise 
            error('No valid primary was selected');
    end
    
    hd = parking_orbit.Altitude;                %Parking orbit altitude
    thetad = parking_orbit.Theta;               %True anomaly at the TTI
    
    %Integrate the stable manifold backwards and check it intersects the whereabouts of the parking orbit
    manifold = 'S';                                                                 %Integrate the stable manifold
    seed = target_orbit.Trajectory;                                                 %Periodic orbit seed
    tspan = target_orbit.tspan;                                                     %Original integration time
    rho = 100;                                                                      %Density of fibres to analyze
    S = invariant_manifold(mu, manifold, branch, seed, rho, tspan, map);            %Initial trajectories
    
    %Relative distance to the primary of interest
    distance = zeros(rho,1);    
    for i = 1:rho
        distance(i) = norm(shiftdim(S.Trajectory(i,S.ArcLength(i),1:3))-R);         %Distance to the primary 
    end
        
    [~, index] = sort(distance);                            %Select the closest manifold to the parking orbit
    s0 = shiftdim(S.Trajectory(index(1),1,:));              %Initial conditions to correct
    Phi = eye(m);                                           %Initial STM 
    Phi = reshape(Phi, [m^2 1]);                            %Initial STM 
    dt = 1e-3;                                              %Time step
    s0 = [s0; Phi];                                         %Initial conditions
    TOF = S.TOF(index(1));                                  %Time of flight
    tspan = TOF:-dt:0;                                                                  %Integration time
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', event);             %Integration tolerances  
    [~, Sn] = ode113(@(t,s)cr3bp_equations(mu, 1, true, t, s), tspan, s0, options);     %Natural trajectory
        
    %Single shooting differential corrector setup 
    GoOn = true;                                                          %Convergence boolean 
    iter = 1;                                                             %Initial iteration
    h_index = 1;                                                          %Altitude range index
    
    if (min(distance) < hd)
        range = linspace(1.05*min(distance), hd, 20);                     %Altitude range to continuate the solution on
    else
        range = linspace(0.95*min(distance), hd, 20);                     %Altitude range to continuate the solution on
    end
    
    %Integration tolerances                  
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', @(t,s)null_flight_path(t, s, R));    %New stopping event
    St = Sn;
    
    %Continuation loop
    while (h_index <= length(range))
        %Preallocation of the maneuver 
        dV = zeros(3,maxIter); 
        
        %Selection of the desired values 
        hd = range(h_index);
        
        %Main correction process
        while ((GoOn) && (iter < maxIter))
            %Compute the error
            Sr = St(end,1:6).'-[R; zeros(3,1)];                             %Relative state vector to the primary of interest
            h = norm(Sr(1:3));                                              %Final altitude
            theta = atan2(Sr(2), Sr(1));                                    %Final 3D true anomaly
                        
            %Error vector
            e = [h-hd];                              
            
            %Dynamics at the end point of the trajectory
            STM = reshape(St(end,m+1:end),[m m]);                           %State transition matrix
            F = cr3bp_equations(mu, 1, false, 0, St(end,1:6).');            %Dynamics vector field
            
            %Compute the sensibility matrix
            dvH = (Sr(1:3)).'/norm(Sr(1:3))*STM(1:3,4:6);                   %Derivative of the altitude with respect to the initial velocity
            dvTheta = ([-Sr(2) Sr(1) 0]/norm(Sr(1:2))^2)*STM(1:3,4:6);      %Derivative of the 3D true anomaly with respect to the initial velocity
            dtH = Sr(1:3).'/norm(Sr(1:3))*F(1:3);                           %Derivative of the altitude with respect to the end time
            dtTheta = ([-Sr(2) Sr(1) 0]/norm(Sr(1:3))^2)*F(1:3);            %Derivative of the 3D true anomaly with respect to the end time
            
            %Derivatives of the flight path angle
            dsGamma = [-Sr(4:6).' -Sr(1:3).'];                              %Derivative with respect to the state
            dtGamma = dsGamma*F(1:6);                                       %Derivative with respect to time
            dsGamma = dsGamma*STM(:,4:6);                                   %Derivative with respect to the velocity
            
            %Sensibility matrix
            Sigma = [dvH]-[dtH]*dsGamma/dtGamma;
            
            %Compute the initial conditions
            dV(:,iter) = -pinv(Sigma)*e;                %Computed maneuver
            s0(4:6) = s0(4:6)+dV(1:3,iter);             %New initial conditions
            
            %Re-integrate the trajectory 
            [~, St] = ode113(@(t,s)cr3bp_equations(mu, 1, true, t, s), tspan, s0, options);
            
            %Convergence analysis 
            if (norm(e) < tol)
                GoOn = false;
            else
                iter = iter+1;
            end          
        end
        
        %Continuation convergence 
        if (iter < maxIter)
            iter = 1;                   %Reinitiate the algorithm
            GoOn = true;                %Reinitiate the algorithm
            h_index = h_index+1;        %Step into the new orbit
        else
            h_index = length(range)+1;  %Stop the continuation process
        end
    end
    
    %Compute the needed maneuver 
    dVf = Sn(end,4:6).'-s0(4:6);                    %Total in-plane velocity change
    
    %Final output 
    xf.Trajectory = St;                             %Final trajectory
    xf.Period = dt*(size(Sn,1)+size(St,1)+1);       %Trajectory time
    state.State = (iter < maxIter);                 %Convergence boolean 
    state.Error = norm(e);                          %Final error
end