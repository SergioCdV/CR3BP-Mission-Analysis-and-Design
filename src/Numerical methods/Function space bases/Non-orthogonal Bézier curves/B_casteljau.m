%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/09/22
% File: B_casteljau.m 
% Issue: 0 
% Validated:

%% B_casteljau %%
% Function for calculation of a Bezier curve using De Casteljau's algorithm 

% Inputs: - array points, the Bézier curve control points 
%         - vector tvec, the control parameter t vector 

% Output: - the resulting Bézier curve

function [curve] = B_casteljau(points,tvec)
    % Number of points
    L = length(points);
    
    % Find number of steps (time increments)
    steps = length(tvec);
    
    % Initialize variable for n-order curve
    curve = zeros(2,steps);
    
    % Use De Casteljau's algorithm to calculate the curve at each interval
    for i = 1:steps
        t = tvec(i);
        P = casteljau(points,L,t);
        curve(:,i) = P(:,1,L);
    end
end