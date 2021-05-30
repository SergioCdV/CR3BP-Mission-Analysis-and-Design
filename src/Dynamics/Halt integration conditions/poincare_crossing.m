%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/05/21
% File: poincare_crossing.m 
% Issue: 0 
% Validated: 30/05/21 

%% Poincaré map crossing %%
% This function contains a handle function to define orbital events to stop integration.

function [Pos, isterminal, dir] = poincare_crossing(~, x)
    %Event definition
    Pos = x(2);         %X axis crossing
    isterminal = 1; 
    dir = 0;            %Direction of the crossing
end