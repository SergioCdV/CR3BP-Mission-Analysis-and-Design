%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Initial guess function %% 
% Function implementation of the a warming up initial guess if available

function [beta, t0, tf] = InitialGuess(params, initial, final)
    % Constants
    mu = params(1); 
    T = params(2); 

    % Approximation of the time of flight via the Vis-Viva theorem
    E0 = jacobi_constant(mu, initial);                         % Energy at the initial state
    Ef = jacobi_constant(mu, final);                           % Energy at the final state
    dE = Ef-E0;                                                % Change of energy
    rho = norm(final(1:3)-initial(1:3));                       % Mean Keplerian radius
    t0 = 0;                                                    % Initial time of flight
    tf = 2*abs(dE)/(T*rho);                                    % Time of flight

    % Preliminary number of revolutions
    dtheta = final(2)-initial(2);
    if (dtheta < 0)
        dtheta = dtheta + 2*pi; 
    end
    
    Napp = ceil( (dtheta+tf*0.5*(initial(4)+final(4)) ) / (2*pi) );
    if (Napp > 0)
        % New initial TOF
        tf = tf*Napp;
    end 

    t0 = 0;
    tf = 20;
    beta = final(2)+2*pi;
end