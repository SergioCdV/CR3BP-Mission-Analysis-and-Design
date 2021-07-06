%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/07/21
% File: zv_plot.m 
% Issue: 0 
% Validated: 

%% Zero velocity surfaces plot %%
% For a given gravitational parameter mu, this function plots in 2D or 3D the zero velocity curves
% corresponding to a given value of the Jacobi Constant.

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - scalar C, the given Jacobi constant
%         - scalar dim, the dimension of the plot

% Outputs:

% New versions: 

function zv_plot(mu, C, dim) 
    if (dim == 3)
        %3D configuration space meshgrid
        d = -2:0.05:2;
        [X,Y,Z] = meshgrid(d, d, d);

        %Compute the pseudo-potential function
        R1 = sqrt((X+mu).^2+Y.^2+Z.^2);     %Relative position to the first primary
        R2 = sqrt((X-1+mu).^2+Y.^2+Z.^2);   %Relative position to the second primary
        Uc = -(X.^2+Y.^2);                  %Centrifugal potential function
        Ug = -2*((1-mu)./R1 + mu./R2);      %Gravitational potential function
        U = Uc+Ug;                          %Total potential function
        
        %Plot results
        figure
        view(3)
        patch(isosurface(X,Y,Z,U,C),'FaceColor','green','EdgeColor','none');
        camlight; lighting phong
        alpha(0.6)
        axis off
        axis equal
        
    else
        %3D configuration space meshgrid
        d = -2:0.01:2;
        [X,Y] = meshgrid(d, d);

        %Compute the pseudo-potential function
        R1 = sqrt((X+mu).^2+Y.^2);          %Relative position to the first primary
        R2 = sqrt((X-1+mu).^2+Y.^2);        %Relative position to the second primary
        Uc = -(X.^2+Y.^2);                  %Centrifugal potential function
        Ug = -2*((1-mu)./R1 + mu./R2);      %Gravitational potential function
        U = Uc+Ug;                          %Total potential function
        
        %Plot results         
        figure
        map = ones(2,3)*0.8;
        colormap(map);
        hold on
        contourf(X, Y, U, C + [0 -0.0001]);
        scatter(-mu, 0, 'filled');
    	scatter(1-mu, 0, 'filled');
        hold off
        axis equal
        axis square
        xlabel('Synodic $x$ coordinate');
        ylabel('Synodic $y$ coordinate');
        title(sprintf('Zero velocity curves corresponding to C = %.2f', C));
    end
end