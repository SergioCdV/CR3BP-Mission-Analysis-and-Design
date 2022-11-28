%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/05/21
% File: control_effort.m 
% Issue: 0 
% Validated: 25/05/21

%% Control effort %%
% This script contains the function to compute the control effort (Lp norm) of a given control law

% Inputs: - vector tspan, the control time span 
%         - 3xm array  u, the control law to be evaluated
%         - boolean discrete_flag, to account for impulsive control laws

% Output: - vector effort, containing the main Lp norms of the control law

% New versions: 

function [effort] = control_effort(tspan, u, discrete_flag)
    % Preallocation 
    effort = zeros(3,1);                          % Figures of merit 
   
    if (discrete_flag)
        % Control effort in time 
        V = sqrt(dot(u,u,1));                               % Euclidean norm of the impulses

        effort(1) = sum(sqrt(dot(u,u,1)));                  % dV integral of the control
        effort(2) = sum(sum(abs(u),1));                     % L1 integral of the control
        effort(3) = sum(dot(u,u,1));                        % L2 integral of the control

        % Infinity norms 
        effort(4) = min(V(V~=0));                           % Minimum impulse
        effort(5) = max(V(V~=0));                           % Maximum impulse
    else
        % Control effort in time 
        a = sqrt(dot(u,u,1));                               % Euclidean norm of the acceleration

        effort(1) = trapz(tspan, a);                        % dV integral of the control
        effort(2) = trapz(tspan, sum(abs(u),1));            % L1 integral of the control
        effort(3) = trapz(tspan, dot(u,u,1));               % L2 integral of the control

        % Infinity norms 
        effort(4) = min(a(a~=0));                           % Minimum acceleration
        effort(5) = max(a(a~=0));                           % Maximum acceleration
    end
end