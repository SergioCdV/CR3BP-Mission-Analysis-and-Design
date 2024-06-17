%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 17/06/24
% File: LibrationPoints.m 
% Issue: 0 
% Validated: 

%% Libration points %%
% For a given gravitational parameter mu, the libration_points function computes 
% the normalized position in the synodic frame of the five libration points of the 
% CR3BP system, using a Newton method for the collinear Li's

% Inputs: - scalar mu, the reduced gravitational parameter of the system 
%         - array R, the position of the primaries in the synodic frame, as a 3x2 array

% Outputs: - structure Lp, containin the definition of each libration point
%          - array L, containing in one matrix the position of the five libration 
%            points of the system (with each column being a position vector) in L1, L2, 
%            L3, L4 and L5 order, as well as the collinear gamma distance
%            to the second primary

% New versions: 

function [Lp, L] = LibrationPoints(mu, R)
    % Preallocation 
    dn = zeros(1,3);                         % Newton step residual
    
    % Compute the equilateral points (forming two symmetric equilateral triangles with the primaries)
    alpha = pi/3;                            % Angle between the libration points and primaries             
    primR = repmat(R(:,1), 1, 2);            % Position of the first primary

    % Equilateral libration points positions
    equiL = primR + [cos(alpha) cos(alpha); sin(alpha) -sin(alpha); 0 0];
    equiL = [equiL; sqrt( dot(equiL-R(:,2), equiL-R(:,2), 1) )];
    
    % Set up the Newton loop for L1/L2/L3
    rh = mu^(1/3);                                              % Hill radius
    lambda = rh * ( 1 + ((-1).^(1:2)) * ( rh/3 + rh^(2/9) ) );  % Initial guess for L1/L2
    rh = 1 - (7/12) * mu;                                       % Initial guess for L3
    lambda = [lambda rh];                                       % Initial guess for L1/L2/L3

    tol = 1E-15;                                                % Newton method tolerance
    iterMax = 1E2;                                              % Maximum allowed iterations for the Newton method
    GoOn = true;                                                % Convergence flag
    iter = 1;                                                   % Initial iteration
        
    % Main computation
    idx = 1:2; 

    while ( GoOn && (iter < iterMax) )
        % Newton algorithm
        f(1:2) = lambda(idx).^5 + (3-mu) * ((-1).^idx) .* lambda(idx).^4 + (3-2*mu) * lambda(idx).^3 - mu * lambda(idx).^2 + 2 *  mu * ((-1).^(idx+1)) .* lambda(idx) - mu;    
        df(1:2) = 5 * lambda(idx).^4 + 4 * (3-mu) * ((-1).^idx) .* lambda(idx).^3 + 3 * (3-2*mu) * lambda(idx).^2 - 2 * mu * lambda(idx) + 2* ((-1).^(idx+1)) * mu;

        f(3) = lambda(3)^5 + (2+mu) * lambda(3)^4 + (1+2*mu) * lambda(3)^3 - (1-mu) * lambda(3)^2 - 2 * (1-mu) * lambda(3) - (1-mu);
        df(3) = 5 * lambda(3)^4 + 4 * (2+mu) * lambda(3)^3 + 3 * (1+2*mu) * lambda(3)^2 - 2 * (1-mu) * lambda(3) - 2 * (1-mu);
    
        dn = -f ./ df;     % Newton step

        % Newton update
        lambda = lambda + dn;
        
        % Check for convergence
        if ( max( abs(dn) ) < tol )
            GoOn = false;
        else
            iter = iter+1; 
        end
    end
    
    % Save the converged collinear point position in an array
    colL = [(1-mu) + ((-1).^(1:2)) .* lambda(1:2) -(mu+lambda(end)); zeros(2,3); lambda];
            
    % Save results in the ouput 
    L = [colL equiL];

    % Create the structure 
    Lp.ID = 1:5;                    % Index of each libration point
    Lp.r = L(1:3,:);                % Position vector of all libration points in Howell's synodic frame
    Lp.gamma = L(end,:);            % Distance to the second primary in normalized units
    Lp.tol = [dn zeros(1,2)];
end