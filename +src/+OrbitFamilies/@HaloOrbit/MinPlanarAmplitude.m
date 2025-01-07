%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 07/01/25
% File: MinPlanarAmplitude.m 
% Issue: 0 

%% Minimum planar amplitude %% 
% This function estimates the minimum halo planar amplitude of a given
% system

% Inputs: - double mu, the reduced gravitational parameter of the system
%         - double libration point L, 1 for L1, 2 for L2 and 3 for L3
%         - double gamma, distance from the libration point to the least
%           massive primary

% Output: - double Ax, the minimum planar amplitude of the system 

function [Ax] = MinPlanarAmplitude(mu, L, gamma)
    % Legendre polynomial coefficients c_n for the Richardson approximation
    order = 4;                                                                  % Order of the approximation
    cn = src.Systems.CR3BPSystem.LegendreCoefficients(mu, L, gamma, order);     % Legendre coefficients
    
    % Determine the orbit spatial eigenvalue    
    polylambda = [1 0 (cn(3) - 2) 0 -(cn(3) - 1) * (1 + 2 * cn(3))];
    lambda = roots( polylambda );

    if ( L == 3 )
        lambda = abs( lambda(3) ) ;
    else        
        lambda = abs( lambda(1) ) ;
    end

    % Richardson 3rd order approximation coefficients
    k = 2 * lambda / (lambda^2 + 1 - cn(3));
    
    d1 = (3 * lambda^2 / k) * (k * (6 * lambda^2 - 1) - 2 * lambda);
    a21 = 3 * cn(4) * (k^2 - 2) / ( 4 * (1 + 2 * cn(3)) );
    a23 = -(3 * cn(4) * lambda / (4 * k * d1)) * (3 * k^3 * lambda - 6 * k * (k - lambda) + 4);
    b21 = -3 * cn(4) * lambda / (2 * d1) * (3 * k * lambda - 4);
    d21 = -cn(4) / (2 * lambda^2);
    s1 = (1.5 * cn(4) * (2 * a21 * (k^2 - 2) - a23 * (k^2 + 2) - 2 * k * b21) - 0.375 * cn(5) * (3 * k^4 - 8 * k^2 + 8)) / (2 * lambda * (lambda * (1 + k^2) - 2 * k));
    a1 = -1.5 * cn(4) * (2 * a21 + a23 + 5 * d21) - 0.375 * cn(5) * (12 - k^2);
    l1 = a1 + 2 * lambda^2 * s1;

    Delta = lambda^2 - cn(3);
    Ax = sqrt( abs( Delta / l1 ) ) * gamma;
end