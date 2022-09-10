%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/05/21
% File: HCNC_guidance.m 
% Issue: 0 
% Validated: 30/05/21

%% Homoclinic Rendezvous Guidance %%
% This script contains the function to compute the guidance law by means of
% an homoclinic connection

% Inputs: - scalar mu, the gravitational parameter of the system
%         - vector branch, selecting the manifold branches to propagate
%         - scalar rho, the number of manifold fiber to compute
%         - scalar L, the libration point ID of the periodic orbit
%         - array target_orbit, with the needed periodic orbit on which the
%           connection is desired to exist
%         - scalar TOF, the time of flight to achieve the rendezvous
%         - boolean long_rendezvous, to allow for primary connections
%         - boolean position fixed, in case the final state is constrained
%         - boolean graphics, to plot or not to plot the results

% Output: - array Sg, the homoclinic trajectory 
%         - scalar dV, the needed maneuver change

% New versions: 

function [Sg, dV] = HCNC_guidance(mu, Branch, rho, L, target_orbit, TOF, long_rendezvous, position_fixed, graphics)
    %General parameters 
    if (position_fixed)
        sd = target_orbit.TargetState;             %Target rendezvous state
    end

    dt = 1e-3;                                     %Time step
    tspan = 0:dt:TOF;                              %Integration time
    target_orbit = target_orbit.Trajectory;        %Target orbit

    %Branch between long/short time rendezvous 
    if (long_rendezvous)  
        %Homoclinic rendezvous Poincaré map
        if (L == 1)
            map = 'Right homoclinic connection';
        else
            map = 'Left homoclinic connection';
        end
        
        %Globalize the unstable manifold until the Poincare section
        manifold = 'U'; 
        branch = Branch(1); 
        Mu = invariant_manifold(mu, L, manifold, branch, target_orbit, rho, tspan, map);
        
        %Globalize the stable manifold until the Poincare section
        manifold = 'S'; 
        branch = Branch(2); 
        Ms = invariant_manifold(mu, L, manifold, branch, target_orbit, rho, tspan, map);

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
            view(3)
            hold on
            plot3(target_orbit(:,1), target_orbit(:,2), target_orbit(:,3), 'b')
            for i = 1:size(Mu.Trajectory,1)
                u = plot3(shiftdim(Mu.Trajectory(i,1:Mu.ArcLength(i),1)),shiftdim(Mu.Trajectory(i,1:Mu.ArcLength(i),2)),...
                      shiftdim(Mu.Trajectory(i,1:Mu.ArcLength(i),3)), 'r');
                u.Color(4) = 0.2;
            end
            for i = 1:size(Ms.Trajectory,1)
                s = plot3(shiftdim(Ms.Trajectory(i,1:Ms.ArcLength(i),1)),shiftdim(Ms.Trajectory(i,1:Ms.ArcLength(i),2)),...
                      shiftdim(Ms.Trajectory(i,1:Ms.ArcLength(i),3)), 'g');
                s.Color(4) = 0.2;
            end
            hold off
            grid on; 
            xlabel('$x$')
            ylabel('$y$')
            zlabel('$z$')
        end

        %Plot the map and select the initial conditions 
        [S0, TOF] = poincare_map(Mu, Ms, 'Connection');
        
        %Integration setup 
        switch (map)
        case 'Left homoclinic connection' 
            map = @(t,s)homoclinic_crossing(t, s, mu, 0);
        case 'Right homoclinic connection'
            map = @(t,s)homoclinic_crossing(t, s, mu, 0);
        end
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', map);
        
        %Integrate the complete trajectory 
        tspan = 0:dt:2*TOF(1);                %Unstable branch time of flight
        
        [tu, Sgu] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, S0.Unstable, options);
        
        %Integration setup         
        tspan = TOF(2):-dt:0;                 %Stable branch time of flight
        [ts, Sgs] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, S0.Stable, options);

        %Reverse the trajectory on the stable manifold
        Sgs = flip(Sgs,1);
        ts = flip(ts);
        
    else
        %Homoclinic rendezvous Poincaré map near the orbit
        map = 'X crossing';            
        
        %Globalize the unstable manifold until the Poincare map
        manifold = 'U'; 
        branch = Branch(1); 
        Mu = invariant_manifold(mu, L, manifold, branch, target_orbit, rho, tspan, map);
        
        %Globalize the stable manifold until the Poincare map
        manifold = 'S'; 
        branch = Branch(2); 
        Ms = invariant_manifold(mu, L, manifold, branch, target_orbit, rho, tspan, map);

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
            view(3)
            hold on
            plot3(target_orbit(:,1), target_orbit(:,2), target_orbit(:,3), 'b')
            for i = 1:size(Mu.Trajectory,1)
                u = plot3(shiftdim(Mu.Trajectory(i,1:Mu.ArcLength(i),1)),shiftdim(Mu.Trajectory(i,1:Mu.ArcLength(i),2)),...
                      shiftdim(Mu.Trajectory(i,1:Mu.ArcLength(i),3)), 'r');
                u.Color(4) = 0.2;
            end
            for i = 1:size(Ms.Trajectory,1)
                s = plot3(shiftdim(Ms.Trajectory(i,1:Ms.ArcLength(i),1)),shiftdim(Ms.Trajectory(i,1:Ms.ArcLength(i),2)),...
                      shiftdim(Ms.Trajectory(i,1:Ms.ArcLength(i),3)), 'g');
                s.Color(4) = 0.2;
            end
            hold off
            grid on; 
            title('Globalization of the unstable and stable manifolds');
            xlabel('Synodic $x$ coordinate')
            ylabel('Synodic $y$ coordinate')
            zlabel('Synodic $z$ coordinate')
            legend('Periodic orbit', 'Unstable manifold', 'Stable manifold', 'Location', 'northeast');
        end
        
        %Plot the map and select the initial conditions 
        [S0, TOF] = poincare_map(Mu, Ms, map);
        
        %Integration setup 
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', @(t,s)poincare_crossing(t,s,mu));
        
        %Integrate the complete trajectory 
        tspan = 0:dt:2*TOF(1);                %Unstable branch time of flight
        
        [tu, Sgu] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, S0.Unstable, options);
        
        %Integration setup         
        tspan = TOF(2):-dt:0;                 %Stable branch time of flight
        [ts, Sgs] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, S0.Stable, options);

        %Reverse the trajectory on the stable manifold
        Sgs = flip(Sgs,1);
        ts = flip(ts);
    end

    %Output 
    Sg.Trajectory = [Sgu; Sgs];                         %Save the trajectory
    Sg.Time = [tu; tu(end)*ones(length(ts),1)+ts];      %Save the time of flight
    Sg.Index = size(Sgu,1);                             %Save the point at which the manifold changes
    dV = S0.Stable(4:6)-Sgu(end,4:6).';                 %Velocity maneuver
end