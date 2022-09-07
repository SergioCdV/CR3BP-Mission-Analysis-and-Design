%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/05/21
% File: figures_merit.m 
% Issue: 0 
% Validated: 25/05/21

%% Figures of merit %%
% This script contains the function to compute the figures of merit related
% to error tracking performance

% Inputs: - vector tspan, the integration time
%         - array S, the trajectory to evaluate

% Output: - vector merit, the figures of merit 

% New versions: 

function [error, merit] = figures_merit(tspan, S)
    % Preallocation 
    merit = zeros(4,1);                             % Figures of merit evaluating the error performance

    % Dimension sanity check 
    if (size(S,2) < size(S,1))
        S = S.';
    end
    
    % Error in time 
    error = sqrt(dot(S,S,1));

    % Compute the error figures of merit 
    merit(1,1) = trapz(tspan, error.^2);              % Integral of the Squared Error                             
    merit(2,1) = trapz(tspan, abs(error));            % Integral of the Absolute value of the Error
    merit(3,1) = trapz(tspan, tspan.*error.^2);       % Integral of the time multiplied by the square error    
    merit(4,1) = trapz(tspan, tspan.*abs(error));     % Integral of the time multiplied by the absolute value of the error
end