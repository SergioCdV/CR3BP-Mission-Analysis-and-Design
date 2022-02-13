%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 11/02/2022
% File: libration_change.m 
% Issue: 0 
% Validated: 

%% Libration change %%
% This function contains a handle function to define orbital events to stop integration.

function [Pos, isterminal, dir] = libration_change(~, mu, Lr, x)
    %Event definition
    L = libration_points(mu);                                   %System libration points
    for i = 1:5
        C(i) = libration_potential(mu,[i;L(:,i)],x,3);
    end
    [~,index] = sort(C); 
    Pos = 2-index(1)
    isterminal = 1; 
    dir = 0;             %Direction of the crossing
end