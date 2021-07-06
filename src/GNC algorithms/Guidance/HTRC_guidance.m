%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/05/21
% File: HTRC_guidance.m 
% Issue: 0 
% Validated: 30/05/21

%% Heteroclinic Rendezvous Guidance %%
% This script contains the function to compute the guidance law by means of
% an heteroclinic connection between periodic orbits

% Inputs: - scalar mu, the gravitational parameter of the system
%         - vector sequence, indicating from which orbit to which orbit the
%           connection is desired
%         - scalar rho, the number of manifold fiber to compute
%         - array target_orbit, the final periodic orbit of the connection
%         - array initial_orbit, the initial periodic orbit of the
%           connection
%         - scalar TOF, the time of flight to achieve the rendezvous
%         - boolean position fixed, in case the final state is constrained
%         - boolean graphics, to plot or not to plot the results

% Output: - array Sg, the homoclinic trajectory 
%         - scalar dV, the needed maneuver change

% New versions: 

function [Sg, dV] = HTRC_guidance(mu, sequence, rho, target_orbit, initial_orbit, TOF, position_fixed, graphics)
    %General parameters 
    if (position_fixed)
        sd = target_orbit.TargetState;             %Target rendezvous state
    end

    dt = 1e-3;                                     %Time step
    tspan = 0:dt:TOF;                              %Integration time
    target_orbit = target_orbit.Trajectory;        %Target orbit
    initial_orbit = initial_orbit.Trajectory;      %Initial orbit

    %Globalize the unstable manifold until the Poincare map
    manifold = 'U'; 
    if (sequence(1) == 1)
        branch = 'R';
        mapu = 'Left heteroclinic connection';
    else
        branch = 'R';
        mapu = 'Right heteroclinic connection';  
    end
    Mu = invariant_manifold(mu, manifold, branch, initial_orbit, rho, tspan, mapu);

    %Globalize the stable manifold until the Poincare map
    manifold = 'S'; 
    if (sequence(2) == 1)
        branch = 'L';
        maps = 'Left heteroclinic connection';
    else
        branch = 'L';
        maps = 'Right heteroclinic connection'; 
    end
    Ms = invariant_manifold(mu, manifold, branch, target_orbit, rho, tspan, maps);

    %Restrict the final position if desired
    if (position_fixed)
        GoOn = true;
        i = 1;
        error = zeros(size(Ms.Trajectory,1),1);
        error(1) = 1e5;
        while (GoOn)
            %Look for the closest fiber to the desired rendezvous position
            error(i+1) = norm(shiftdim(Ms.Trajectory(i,1,:))-sd);
            if (error(i+1) < error(i))
                Trajectory = Ms.Trajectory(i,:,:); 
                ArcLength = Ms.ArcLength(i);
                TOF = Ms.TOF(i);
                Index = Ms.Index(i);
                GoOn = false;
            else
                if (i+1 > size(Ms.Trajectory,1))
                    GoOn = false;
                else
                    i = i+1;
                end
            end
        end

        %Save results
        Ms.Trajectory = Trajectory;
        Ms.ArcLength = ArcLength; 
        Ms.TOF = TOF;
        Ms.Index = Index;
    end

    %Plot the manifold
    if (graphics)
        figure
        hold on
        plot3(target_orbit(:,1), target_orbit(:,2), target_orbit(:,3), 'b')
        plot3(initial_orbit(:,1), initial_orbit(:,2), initial_orbit(:,3), 'k')
        for i = 1:size(Mu.Trajectory,1)
            plot3(shiftdim(Mu.Trajectory(i,1:Mu.ArcLength(i),1)),shiftdim(Mu.Trajectory(i,1:Mu.ArcLength(i),2)),...
                  shiftdim(Mu.Trajectory(i,1:Mu.ArcLength(i),3)), 'r');
        end
        for i = 1:size(Ms.Trajectory,1)
            plot3(shiftdim(Ms.Trajectory(i,1:Ms.ArcLength(i),1)),shiftdim(Ms.Trajectory(i,1:Ms.ArcLength(i),2)),...
                  shiftdim(Ms.Trajectory(i,1:Ms.ArcLength(i),3)), 'g');
        end
        hold off
        grid on; 
        title('Heteroclinic connections using invariant manifolds');
        xlabel('Synodic $x$ coordinate')
        ylabel('Synodic $y$ coordinate')
        zlabel('Synodic $z$ coordinate')
        legend('Target orbit', 'Initial orbit', 'Unstable manifold', 'Stable manifold', 'Location', 'northeast');
    end

    %Plot the map and select the initial conditions 
    [S0, TOF] = poincare_map(Mu, Ms, 'Connection');

    %Integration setup 
    switch (mapu)
        case 'Left heteroclinic connection' 
            map = @(t,s)heteroclinic_crossing(t, s, mu, 1);
        case 'Right heteroclinic connection'
            map = @(t,s)heteroclinic_crossing(t, s, mu, -1);
    end
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', map);

    %Integrate the complete trajectory 
    tspan = 0:dt:2*TOF(1);                %Unstable branch time of flight

    [tu, Sgu] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, S0.Unstable, options);

    %Integration setup     
    switch (maps)
    case 'Left heteroclinic connection' 
        map = @(t,s)heteroclinic_crossing(t, s, mu, 1);
    case 'Right heteroclinic connection'
        map = @(t,s)heteroclinic_crossing(t, s, mu, -1);
    end
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', map);
    tspan = TOF(2):-dt:0;                 %Stable branch time of flight
    
    [ts, Sgs] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, S0.Stable, options);

    %Reverse the trajectory on the stable manifold
    Sgs = flip(Sgs,1);
    ts = flip(ts);

    %Output 
    Sg.Trajectory = [Sgu; Sgs];                         %Save the trajectory
    Sg.Time = [tu; tu(end)*ones(length(ts),1)+ts];      %Save the time of flight
    Sg.Index = size(Sgu,1);                             %Save the point at which the manifold changes
    dV = S0.Stable(4:6)-Sgu(end,4:6).';                 %Velocity maneuver
end
