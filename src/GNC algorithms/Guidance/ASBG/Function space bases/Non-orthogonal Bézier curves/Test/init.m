%% Project: 
% Date: 29/01/22

%% Init
% This function initializes graph for Bezier curve. Plots given points and lines connectiong them. Option to increase graph margins

function init(points, margin)
    % Find edge points for axis scaling
    M = max(points, [], 2);
    m = min(points, [], 2);
    
    % Define axes and margins
    axis([(m(1)-margin) (M(1)+margin) (m(2)-margin) (M(2)+margin)]);
    
    % Colour for polygon lines
    c = [210, 211, 214]./255;
    
    % Plot the initial points and their connecting lines
    plot(points(1,:),points(2,:),'Color',c,'LineWidth',3);
    scatter(points(1,:),points(2,:),'filled','k');
    hold on
end