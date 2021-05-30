%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/05/21
% File: poincare_map.m 
% Issue: 0 
% Validated: 30/05/21

%% Poincare map %%
% This script contains the function to compute several Poincare surfaces of
% section

% Inputs: - structure safe_corridor, the parameters of the safety policy 
%         - structure Penalties, the controller parameters
%         - array So, the relative position of the obstacles 
%         - scalar TOF, the time of flight to achieve the rendezvous
%         - vector s0, the initial conditions of the relative particle
%         - string map, selecting the Poincaré map to plot

% Output: - array Sg, the guidance law definition parameters

% New versions: 

function [S0, dV, TOF] = poincare_map(Su, Ss, map)
    %Branch the different Poincarñe maps
    switch (map)
        case 'Homoclinic rendezvous'
            %Extract the intersection points 
            Iu = shiftdim(Su.Trajectory(:,Su.ArcLength,:));    %Intersection of the unstable manifold
            Is = shiftdim(Ss.Trajectory(:,Ss.ArcLength,:));    %Intersection of the stable manifold
            
            %Select the connection point in the configuration space
            figure(1) 
            hold on 
            for i = 1:size(Iu,1)
                plot(Iu(i,1,1), Iu(i,1,3), 'ro');
            end
            for i = 1:size(Is,1)
                plot(Is(i,1,1), Is(i,1,3), 'go');
            end
            hold off
            grid on; 
            %[epoch, RAANs] = getpts();   
            
            %Extract the initial conditions on the unstable manifold
            
            %Determine the size of the impulsive maneuver
            
        otherwise 
            error('No valid map was selected');
    end
end