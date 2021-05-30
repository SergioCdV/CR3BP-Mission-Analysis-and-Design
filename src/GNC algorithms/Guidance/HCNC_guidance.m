%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/05/21
% File: HCNC_guidance.m 
% Issue: 0 
% Validated: 30/05/21

%% Homoclinic Rendezvous Guidance %%
% This script contains the function to compute the guidance law by means of
% an homoclinic connection

% Inputs: - structure safe_corridor, the parameters of the safety policy 
%         - structure Penalties, the controller parameters
%         - array So, the relative position of the obstacles 
%         - scalar TOF, the time of flight to achieve the rendezvous
%         - vector s0, the initial conditions of the relative particle
%         - boolean offline_flag, to allow an offline computation of the
%           guidance law

% Output: - array Sg, the guidance law definition parameters

% New versions: 

function [Sg, dV] = HCNC_guidance(mu, Branch, rho, seed, TOF, long_rendezvous)
    %Branch between long/short time rendezvous 
    if (long_rendezvous)
        %General parameters 
        dt = 1e-3;                     %Time step
        tspan = 0:dt:seed.Period;      %Integration time
        seed = seed.Trajectory;        %Target orbit
        
        map = 'Homoclinic rendezvous';         %Homoclinic rendezvous Poincaré map
        
        %Globalize the unstable manifold until the Poincare map
        manifold = 'U'; 
        branch = Branch(1); 
        Mu = invariant_manifold(mu, manifold, branch, seed, rho, tspan, map);
        
        %Globalize the stable manifold until the Poincare map
        manifold = 'S'; 
        branch = Branch(2); 
        Ms = invariant_manifold(mu, manifold, branch, seed, rho, tspan, map);
        
        %Plot the map and select the initial conditions 
        [s0, dV, TOF] = poincare_map(Mu, Ms, map);
        
        %Integration setup 
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', @(t,s)poincare_crossing(t,s));
        
        %Integrate the complete trajectory 
        tspan = 0:dt:TOF(1);                %Unstable branch time of flight
        
        [tu, Sgu] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, s0, options);
        
        %Integration setup 
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
        
        tspan = 0:dt:TOF(2);                %Stable branch time of flight
        s0 = Sgu(end,:);                    %Stable manifold initial conditions
        s0(4:6) = s0(4:6) + dV;             %Stable manifold initial conditions
        
        [ts, Sgs] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, s0, options);
        
        %Output 
        Sg.Trajecory = [Sgu(1:end-1,:); Sgs];
        Sg.Time = [tu; tu(end)*ones(length(ts),1)+ts];
    else
    end
end