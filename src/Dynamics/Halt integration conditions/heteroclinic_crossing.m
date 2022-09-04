%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/06/21
% File: heteroclinic_crossing.m 
% Issue: 0 
% Validated: 06/07/21 

%% Poincar� map crossing %%
% This function contains a handle function to define orbital events to stop integration

function [Pos, isterminal, dir] = heteroclinic_crossing(~, x, mu, sign)
    % Event definition
    Pos = x(1)-(1-mu);      % Crossing of the Poincar� map at the secondary primary
    isterminal = 1;         % Halt the integration
    dir = sign;             % Direction of the crossing
end