%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: libration_points.m 
% Issue: 0 
% Validated: 

%% Libration points %%
% For a given gravitational parameter mu, the libration_points function
% computes the normalized position in the synodic frame of the five libration points of the CR3BP system,
% using a Newton method for the collinear Li's.

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 

% Outputs: - array L, containing in one matrix the position of the five libration 
%            points of the system (with each column being a position vector) in L1, L2, 
%            L3, L4 and L5 order, as well as the collinear gamma distance for the first three.

% New versions: 

function [L] = libration_points(mu)
    %Preallocation 
    colL = zeros(4,3); 
    
    %Compute the equilateral points (forming two symmetric equilateral triangles with the primaries)
    alpha = pi/3;                               %Angle between the libration points and primaries             
    primR = [-mu -mu; 0 0; 0 0; 1 1];           %Position of the first primary
    equiL = primR+[cos(alpha) cos(alpha); 
                   sin(alpha) -sin(alpha); 
                                  0 0; 0 0];    %Libration points positions
    
    %Compute the collinear points using a standard Newton method
    numL = 2;           %Number of collinear points to calculate recursively
    tol = 1e-15;        %Newton method tolerance
    iterMax = 1e5;     %Maximum allowed iterations for the Newton method
    
    %Main loop
    for i = 1:numL
        %Set up the Newton loop for L1/L2
        rh = mu^(1/3);                              %Hill radius
        lambda = rh*(1+((-1)^i)*(rh/3 +rh^2/9));    %Initial guess 
        GoOn = true;                                %Convergence flag
        iter = 1;                                   %Initial iteration
        
        %Main computation
        while ((GoOn) && (iter < iterMax))
            %Newton algorithm
            f = lambda^5 +((-1)^i)*(3-mu)*lambda^4 +(3-2*mu)*lambda^3 ...
                -mu*lambda^2 +2*((-1)^(i+1))*mu*lambda-mu;
            df = 5*lambda^4 +((-1)^i)*4*(3-mu)*lambda^3 +3*(3-2*mu)*lambda^2 ...
                 -2*mu*lambda +2*((-1)^(i+1))*mu;
            dn = -f/df;
            lambda = lambda +dn;
            
            %Check for convergence
            if (abs(dn) < tol)
                GoOn = false;
            else
                iter = iter+1; 
            end
        end
        
        %Save the converged collinear point position in an array
        colL(:,i) = [(1-mu)+((-1)^(i))*lambda(end); 0; 0; lambda];
    end
    
    %Set up the Newton loop for L3
    rh = 1-(7/12)*mu;       %Initial guess
    lambda = rh;            %Initial guess
    GoOn = true;            %Convergence flag
    iter = 1;               %Initial iteration
    
    %Main computation
    while ((GoOn) && (iter < iterMax))
        f = lambda^5 +(2+mu)*lambda^4 +(1+2*mu)*lambda^3 -(1-mu)*lambda^2 -2*(1-mu)*lambda -(1-mu);
        df = 5*lambda^4 +4*(2+mu)*lambda^3 +3*(1+2*mu)*lambda^2 -2*(1-mu)*lambda -2*(1-mu); 
        dn = -f/df;
        lambda = lambda +dn;
        if (abs(dn) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end
    colL(:,3) = [-(mu+lambda); 0; 0; lambda];
    
    %Save results in the L structure ouput 
    L = [colL equiL];
end