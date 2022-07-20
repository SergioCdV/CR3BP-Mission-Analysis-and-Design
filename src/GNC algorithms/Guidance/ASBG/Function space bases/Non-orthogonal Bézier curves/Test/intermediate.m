%% Project: 
% Date: 29/01/22

%% Intermediate
% This function uses De Casteljau's algorithm to  evaluate the set of points in each instant "t", and calculate the
% auxilairy points and lines. The outputs to the function are the plotted data, which can be deleted
% to clear graph.

function [linemarker, dotmarker] = intermediate(points,t)
    % Number of points
    L = size(points,2);
    
    % Use time variable as marker
    xlabel(sprintf('t = %.2f', t));
    
    % Generate colour gradient for plots
    c = parula(L);
    
    % Calulate Bezier curve (recursive)
    P = casteljau(points,L,t);
    
    %Preallocation for speed 
    dotmarker = zeros(L,L); 
    linemarker = zeros(L,L); 
    
    % Plotting in a double iteration
    for i = 1:L
        for j = 1:L-i     
            % Plot auxiliary points
            if i==(L-1)
                dotmarker(i,j) = scatter(P(1,j,i+1),P(2,j,i+1),40,c(i,:),'filled');
            else
                dotmarker(i,j) = scatter(P(1,j,i+1),P(2,j,i+1),10,c(i,:));
            end
    
            % Plot the lines between points
            if (j > 1) && (L > 2)
                linemarker(i,j) = plot(P(1,j-1:j,i+1),P(2,j-1:j,i+1),'Color',c(i,:));
            elseif (L < 3)
                linemarker = [];
            end
        end
    end
end