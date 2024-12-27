%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 07/05/21
% File: fp_crossing.m 
% Issue: 0 
% Validated: 

%% First primary crossing %%
% This function contains a handle function to define orbital events to stop integration

function [Pos, isterminal, dir] = fp_crossing(~, x, mu)
    %Event definition
    Pos = x(1)+mu;                  % Crossing the first primary Poincaré map
    isterminal = 1;                 % Halt the integration
    dir = -1;                       % Direction of the crossing
end