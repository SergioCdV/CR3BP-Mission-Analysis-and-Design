%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/07/21
% File: potential_plot.m 
% Issue: 0 
% Validated: 

%% Potential function plot %%
% For a given gravitational parameter mu, this function plots the potential function as well as 
% the location of its saddle points.

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - scalar dim, the dimension of the plot (2D or 3D)


% Outputs:

% New versions: 

function potential_plot(mu, dim) 
    %Planar configuration space meshgrid 
    x = -2:1e-2:2;                  %Synodic x coordinate
    y = -2:1e-2:2;                  %Synodic y coordinate
    [X,Y] = meshgrid(x,y);          %Complete meshgrid
    
    %Compute the associated potential function in vectorized form
    R1 = sqrt((X+mu).^2+Y.^2);      %Relative position to the first primary
    R2 = sqrt((X-1+mu).^2+Y.^2);    %Relative position to the second primary
    Uc = -(X.^2+Y.^2);              %Centrifugal potential function
    Ug = -2*((1-mu)./R1+mu./R2);    %Gravitational potential function
    U = Uc+Ug;                      %Total potential function
   
    %Plot results 
    figure
    hold on
    if (dim == 3)
        surf(X, Y, U, 'FaceColor', [1 0.99 1], 'EdgeColor', 'none');
        camlight left; lighting phong
        zlim([-3.8 -2.83])
        alpha(0.5)
        view([60 80])
        zlabel('$\tilde{U}$');
    else
        rho = 500;                             %Number of isocurves
        contour(X, Y, U, rho);                 %Isoplot
    end
    scatter(-mu, 0, 'filled');
    scatter(1-mu, 0, 'filled');
    hold off
    grid on;
    xlabel('Synodic $x$ coordinate');
    ylabel('Synodic $y$ coordinate');
    title('Pseudo-potential function $\tilde{U}$ surface')
end