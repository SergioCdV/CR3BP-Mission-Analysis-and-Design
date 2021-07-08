%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 07/07/21
% File: energy_cases.m 
% Issue: 0 
% Validated: 

%% Energy cases plot %%
% For a given gravitational parameter mu, this function plots the prohibited realms associated to 
% the libration points Jacobi constant level.

% Inputs: - scalar mu, the reduced gravitational parameter of the system

% Outputs:

% New versions: 

function energy_cases(mu)
    %Compute the location of the libration points
    L = libration_points(mu); 
    
    %Compute the libration point Jacobi constant levels 
    J = JacobiHill_values(mu);
    
    %Plot the associated isocontours of the pseudo-potential function 
    for i = 1:length(J)-1
        %Plot the realms
        figure(i)
        zv_plot(mu, J(i), 2);
        
        %Mark the libration points
        figure(i) 
        hold on 
        if (i ~= length(J)-1)
            plot(L(1,i), L(2,i), 'ko');
        else
            plot(L(1,i), L(2,i), 'ko');
            plot(L(1,i+1), L(2,i+1), 'ko');
        end
        hold off
    end
end