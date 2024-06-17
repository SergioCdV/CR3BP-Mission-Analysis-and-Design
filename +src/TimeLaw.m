%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 17/06/24
% File: TimeLaw.m 
% Issue: 0 
% Validated: 

%% Time Law %%
% For a given CR3BP, this function computes the motion of the second
% primary with respect to the first in its circular motion

% Inputs: - scalar T, the period of the system
%         - scalar t, the current epoch 
%         - scalar t0, the reference epoch, at which theta is 0

% Outputs: - scalar theta, the angular coordinate of the second primary
%            along the orbit of the first primary

% New versions: 

function [theta] = TimeLaw(T, t0, t)
    % Compute the relative span
    dt = t - t0;                    % Relative epoch
    dt = mod(dt, T);                % Non-dimensional span

    % Output 
    theta = dt;
end