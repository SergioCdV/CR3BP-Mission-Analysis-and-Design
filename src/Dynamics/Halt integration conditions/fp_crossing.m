%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 07/05/21
% File: fp_crossing.m 
% Issue: 0 
% Validated: 

%% First primary crossing %%
% This function contains a handle function to define orbital events to stop integration.

function [Pos, isterminal, dir] = fp_crossing(~, x, mu)
    %Event definition
    Pos = x(1)+mu;                  %First primary crossing
    isterminal = 1; 
    dir = -1;                        %Direction of the crossing
end