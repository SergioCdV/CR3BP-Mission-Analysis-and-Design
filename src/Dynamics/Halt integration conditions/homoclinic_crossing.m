%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/05/21
% File: homoclinic_crossing.m 
% Issue: 0 
% Validated: 30/05/21 

%% Poincaré map crossing %%
% This function contains a handle function to define orbital events to stop integration.

function [Pos, isterminal, dir] = homoclinic_crossing(~, x, mu, sign)
    %Event definition
    if (x(1) > 1-mu)
        Pos = x(2);         %X axis crossing
    else
        Pos = NaN;
    end
    isterminal = 1;         %Halt the integration
    dir = sign;             %Direction of the crossing
end