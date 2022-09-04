%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 25/07/22

%% Sundman radius %%
% Function to compute the Sundman generalized function r

% Inputs: - scalar mu, the gravitational parameter of the system
%         - structure St, defining the target evolution and the vectorfield
%           to be used
%         - vector tau, the regularized sampling grid
%         - array C, the 9xm state vector

% Outputs: - vector r, the Sundman generalized function

function [r] = sundman_radius(mu, sf, St, tau, C)
    % Location of the primaries
    R = repmat([-mu; zeros(5,1); 1-mu; zeros(5,1)],1,length(tau));
    
    % Impose the correct dimensions for the preliminary sampling grid and evaluate the target trajectory if needed
    switch (St.Field)
        case 'Absolute'
            St.Trajectory = repmat(St.Trajectory(:,1),1,length(tau));   
    
        case 'Relative'
            St.Trajectory = target_trajectory(sf, tau, St.Period, St.Cp);
    end
    
    % Relative position vector
    s = cylindrical2cartesian(C, true);                     % Relative Cartesian state vector
    Rr(1:3,:) = s(1:3,:)-R(1:3,:)+St.Trajectory(1:3,:);     % Position to the first moving primary
    Rr(4:6,:) = s(1:3,:)-R(4:6,:)+St.Trajectory(1:3,:);     % Position to the second moving primary 

    % Sundman radius
    r = (1-mu)./sqrt(Rr(1,:).^2+Rr(2,:).^2+Rr(3,:).^2)+mu./sqrt(Rr(4,:).^2+Rr(5,:).^2+Rr(6,:).^2);
end