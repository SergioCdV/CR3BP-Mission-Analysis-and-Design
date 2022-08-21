%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 21/08/22
% File: yrgb_guidance.m 
% Issue: 0 
% Validated: 21/08/22

%% Final orbit %%
% Function to compute the target's insertion condition on the
% quasi-periodic curve

% Inputs: - array curve, the target's orbit QPC Chebyshev polynomial weights
%         - vector theta, the evaluation torus latitudes

% Outputs: - array C, the final target's periodic evolution

function [C] = final_orbit(curve, theta)
    theta = 2*theta/(2*pi)-1;                             % Mapping theta into the 2*pi * [-1,1] domain
    P = CH_basis('first', size(curve,2)-1, theta);        % Chebyshev polynomials
    C = curve*P;                                          % Evaluate the target trajectory
end