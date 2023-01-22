%% Project: 
% Date: 29/01/2022

%% Casteljau
% Function for calculation of Bezier curve using De Casteljau's algorithm

% Inputs: - array points, the Bézier curve control points 
%         - vector tvec, the control parameter t vector 

% Output: - the resulting Bézier curve

% Iterative portion of De Casteljau's algorithm
function [P] = casteljau(points, L, t)
    % Extract points to the first iteration (defined points)
    P(:,:,1) = points;
    
    for i = 1:L
        for j = 1:L-i
            % Calulate Bezier curve (recursively)
            P(:,j,i+1) = (1-t)*P(:,j,i) + t*P(:,j+1,i); 
        end
    end
end