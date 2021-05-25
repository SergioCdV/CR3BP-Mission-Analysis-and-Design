%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/05/21
% File: figures_merit.m 
% Issue: 0 
% Validated: 25/05/21

%% Figures of merit %%
% This script contains the function to compute the figures of merit of the trajectory error.

% Inputs: - vector tspan, the integration time
%         - array S, the trajectory to evaluate

% Output: - vector merit, the figures of merit 

% New versions: 

function [error, merit] = figures_merit(tspan, S)
    %Preallocation 
    merit = zeros(4,1);                             %Figures of merit evaluating the rendezvous trajectory
    error = zeros(length(tspan),1);                 %Error to rendezvous 
    dumb = zeros(length(tspan),1);                  %Just an auxiliary variable
    
    %Error in time 
    for i = 1:length(tspan)
        error(i) = norm(S(i,7:12));
    end

    %Compute the error figures of merit 
    merit(1) = trapz(tspan, error.^2);              %Integral of the square error                             
    merit(2) = trapz(tspan, abs(error));            %Integral of the absolute value of the error
    
    for i = 1:length(tspan)
        dumb(i) = tspan(i)*error(i)^2;
    end
    merit(3) = trapz(tspan, dumb);                  %Integral of the time multiplied by the square error    
    
    for i = 1:length(tspan)
        dumb(i) = tspan(i)*abs(error(i));
    end
    merit(4) = trapz(tspan, dumb);                  %Integral of the time multiplied by the absolute value of the error
end