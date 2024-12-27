%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/2024
% File: CoLegendreCoefficients.m 
% Issue: 0 
% Validated: 

%% Relative Legendre coefficients %%
% This function computes the Legendre coefficients of the Legendre expasion of the relative motion 

% Inputs: - scalar mu, the reduced gravitational parameter of the system 
%         - array r_t [3 x N], the target spacecraft synodic position vector
%         - scalar order, the order of the expansion

% Outputs: - vector cn [order+1xN], containing the Legendre coefficients up to order

function [cn] = CoLegendreCoefficients(mu, r_t, order)
    % Characteristics of the system 
    mup(1) = 1-mu;                         % Reduced gravitational parameter of the first primary
    mup(2) = mu;                           % Reduced gravitational parameter of the second primary
    R(:,1) = [-mu; 0; 0];                  % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];                 % Synodic position of the second primary
    
    % Preallocation of the coefficients 
    Rr(1:3,:) = R(:,1) - r_t;              % Synodic relative position of the target to the first primary
    Rr(4:6,:) = R(:,2) - r_t;              % Synodic relative position of the target to the second primary

    cn = mup(1) ./ sqrt( dot(Rr(1:3,:), Rr(1:3,:), 1) ).^( (0:order).' + 1 ) + mup(2) ./ sqrt( dot(Rr(4:6,:), Rr(4:6,:), 1) ).^( (0:order).' + 1 );
end