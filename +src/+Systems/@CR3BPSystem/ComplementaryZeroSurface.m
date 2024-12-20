%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 19/12/24
% File: ComplementaryZeroSurface.m 
% Issue: 0 
% Validated: 

%% Complementary Zero Velocity Surface %%
% For a given gravitational parameter mu and Jacobi constant C, this
% function computes the region of the configuration space denotes as Complementary Zero
% Velocity Surfaces, as given in Oshima, 2024

% Inputs: - double mu, the reduced gravitational parameter of the system 
%         - double C, the Jacobi constant for which the ZVS is to be
%           computed

% Outputs: - array s [3xN], the regions of the sampled phase space matching
%            the energy constraint

% New versions: 

function [s] = ComplementaryZeroSurface(mu, C, display_flag)
    % Sanity checks 
    if ( ~exist("display_flag", "var") )
        display_flag = false;
    end

    % Constants 
    tol = 1E-3;                                   % Tolerance to the identification of the surface
    R(:,1) = [-mu; 0; 0];                         % Position vector of the first primary
    R(:,2) = [1 - mu; 0; 0];                      % Position vector of the second primary

    % Sample the fundamental XY plane 
    d = linspace(-2, +2, 100);                    % Sampling of the configuration space
    [X, Y, Vz] = meshgrid(d, d, d);               % Sampled configuration space

    % Compute the isocurves 
    R1 = sqrt( (X + mu).^2 + Y.^2 );              % Relative position to the first primary
    R2 = sqrt( (X - (1-mu)).^2 + Y.^2 );          % Relative position to the second primary

    Ug = - (1 - mu) ./ R1 - mu ./ R2;             % Gravitational term
    Uc = - 0.5 * (X.^2 + Y.^2);                   % Centrifugal potential
    U = Ug + Uc;                                  % Total potential
    Jc = -2 * U - Vz.^2;                          % Jacobi Constant
    dC = Jc - C;                                  % Residual

    idx = abs(dC) < tol;
    s = [X(idx).'; Y(idx).'; Vz(idx).'];

    % Plot the results 
    if (display_flag)
        figure
        src.graphics.set_graphics()
        view(3)
        hold on
        patch(isosurface(X, Y, Vz, Jc, C), 'FaceColor', 'red', 'EdgeColor', 'none');
        camlight; lighting phong
        alpha(0.4)
        scatter3(R(1,1), R(2,1), R(3,1), 'k', 'filled');
        scatter3(R(1,2), R(2,2), R(3,2), 'k', 'filled');
        hold off
        labels = {'$M_1$', '$M_2$'};
        text(R(1,:), R(2,:)+0.1, [1e-1 1e-1], labels);
        grid on; 
        xlabel('$x$');
        ylabel('$y$');
        zlabel('$v_z$');
        title(sprintf('CZS @ C = %.2f', C));
    end
end