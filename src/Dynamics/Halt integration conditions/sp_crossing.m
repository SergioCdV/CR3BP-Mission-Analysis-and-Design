%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 07/05/21
% File: sp_crossing.m 
% Issue: 0 
% Validated: 

%% Secondary primary crossing %%
% This function contains a handle function to define orbital events to stop integration

function [Pos, isterminal, dir] = sp_crossing(~, x, mu)    
    %Event definition
    Pos = x(1)-(1-mu);              % Secondary primary crossing
    isterminal = 1;                 % Halt the integration
    dir = 0;                        % Direction of the crossing
end