%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/09/22
% File: B_classic.m 
% Issue: 0 
% Validated:

%% B_classic
% This function uses the classical implementation for calculating a Bezier curve as defined by the determined points

% Inputs: - array points, the Bézier curve control points 
%         - vector tvec, the control parameter t vector 

% Output: - the resulting Bézier curve B

function [B] = B_classic(points,tvec)
    % Determine order
    n = length(points)-1;

    % Generate the polynomial basis 
    bern = bernstein_basis(n,tvec);
    
    % Bézier curve as a dot product
    B = points*bern;
end