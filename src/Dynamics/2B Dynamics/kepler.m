%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/10/21
% File: kepler.m 
% Issue: 0 
% Validated: 12/10/21

%% Kepler %%
% This file contains the function to solve the time law of the Kepler problem, 
% by means of Laguerre-Conway iterations

% Inputs: - vector elements, containing the mean classical Euler orbital elements (a, e, RAAN, i, omega, M) 

% Ouputs: - scalar theta, containing the true anomaly associated to the elements's mean anomaly

% New version updates: 

function [theta] = kepler(elements)
    % Constants 
    e = elements(2);                                            % Eccentricity of the orbit
    M0 = elements(6);                                           % Mean anomaly of the orbit 
    
    % Set up the Newton method 
    k = 5;                                                      % Conway constant
    tol = 1e-15;                                                % Convergence tolerance
    iterMax = 10^6;                                             % Maximum number of iterations
    GoOn = true;                                                % Convergence flag
    iter = 1;                                                   % Initial iteration
    u = M0+e;                                                   % Conway method variable
    E(iter) = (M0*(1-sin(u))+u*sin(M0))/(1+sin(M0)-sin(u));     % Initial guess for the eccentric anomaly
    
    % Main computation 
    while ((GoOn) && (iter < iterMax))
        % Laguerre-Conway iterations
        f = E(iter)-e*sin(E(iter))-M0;                          % Kepler equation
        df = 1-e*cos(E(iter));                                  % First-derivative of the Kepler equation
        ddf = e*sin(E(iter));                                   % Second-derivative of the Kepler equations
        dn = -k*(f)/(df+sqrt((k-1)^2*df^2-k*(k-1)*f*ddf));      % Newton step

        % Newton update
        E(iter+1) = E(iter)+dn;
        
        % Convergence checking 
        if (abs(dn) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end  
    
    % Final tsrue anomaly
    theta = atan2(sqrt(1-e^2)*sin(E(end))/(1-e*cos(E(end))), (cos(E(end))-e)/(1-e*cos(E(end))));
end