%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/05/21
% File: control_effort.m 
% Issue: 0 
% Validated: 25/05/21

%% Control effort %%
% This script contains the function to compute the control effort (Lp norm) of a given control law.

% Inputs: - vector tspan, the integration time 
%         - vector u, the control law to be evaluated 

% Output: - vector effort, containing the main Lp norms of the control law

% New versions: 

function [effort] = control_effort(tspan, u)
    %Preallocation 
    effort = zeros(size(u,1),2);                              %Figures of merit evaluating the rendezvous trajectory
    
    %Error in time 
    for i = 1:size(u,1)
        effort(i,1) = trapz(tspan, u(i,:).^2);                 %L2 integral of the control
        effort(i,2) = trapz(tspan, sum(abs(u(i,:)),1));        %L1 integral of the control
    end
end