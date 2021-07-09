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
    %Set graphics 
    set_graphics(); 
    
    %Planar configuration space meshgrid 
    x = -1.5:1e-2:1.5;              %Synodic x coordinate
    y = -1.5:1e-2:1.5;              %Synodic y coordinate
    [X,Y] = meshgrid(x,y);          %Complete meshgrid
    
    %Compute the associated potential function in vectorized form
    R1 = sqrt((X+mu).^2+Y.^2);      %Relative position to the first primary
    R2 = sqrt((X-1+mu).^2+Y.^2);    %Relative position to the second primary
    Uc = -(X.^2+Y.^2);              %Centrifugal potential function
    Ug = -2*((1-mu)./R1+mu./R2);    %Gravitational potential function
    U = Uc+Ug;                      %Total potential function
    
    %Location of the primaries  
    R(:,1) = [-mu; 0];              %Location of the first primary
    R(:,2) = [1-mu; 0];             %Location of the second primary
   
    %Plot results 
    hold on
    if (dim == 3)
        surf(X, Y, U, 'FaceColor', 'red', 'EdgeColor', 'none');
        camlight headlight; lighting phong
        zlim([-3.5 -2.8])
        view([60 70])
        zlabel('$\tilde{U}$');
        xlabel('Synodic $x$ coordinate'); 
        ylabel('Synodic $y$ coordinate');
        alpha(0.9);
        title('Pseudo-potential function $\tilde{U}$ surface')
        grid on;
    else
        rho = -3.5:5e-2:-2.5;              %Number of isocurves
        contour(X,Y,U,rho);                %Isoplot
        scatter(R(1,1), R(2,1), 'k', 'filled');
        scatter(R(1,2), R(2,2), 'k', 'filled');
        labels = {'$M_1$', '$M_2$'};
        text(R(1,:), R(2,:)+0.1, labels);
        hold off
        grid on;
        colorbar
        xlabel('Synodic $x$ coordinate');
        ylabel('Synodic $y$ coordinate');
        title('Pseudo-potential function $\tilde{U}$ isocurves')
    end
end