%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: x_crossing.m 
% Issue: 0 
% Validated: 

%% x axis crossing %%
% This function contains a handle function to define orbital events to stop integration

function [Pos, isterminal, dir] = x_crossing(~, x)
    % Event definition
    Pos = x(2);         % X axis crossing
    isterminal = 1;     % Halt the integration
    dir = -1;           % Direction of the crossing
end