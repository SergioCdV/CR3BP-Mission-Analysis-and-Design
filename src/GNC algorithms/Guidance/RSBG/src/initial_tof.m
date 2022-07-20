%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 16/04/22

%% Initial TOF %%
% Function to estimate the initial time of flight

% Inputs: - scalar mu, the gravitational parameter of the system  
%         - scalar T, the maximum acceleration allowed for the
%           spacecraft
%         - vector initial, the initial boundary conditions of the
%           trajectory in Cartesian coordinates
%         - vector final, the final boundary conditions of the
%           trajectory in Cartesian coordinates

% Outputs:- scalar tfapp, the initial approximation of the time of flight

function [tfapp] = initial_tof(mu, T, initial, final)
    % Approximation of the time of flight via the Vis-Viva theorem
    E0 = jacobi_constant(mu, initial);                         % Energy at the initial state
    Ef = jacobi_constant(mu, final);                           % Energy at the final state
    dE = Ef-E0;                                                % Change of energy
    rho = norm(final(1:3)-initial(1:3));                       % Mean Keplerian radius
    tfapp = 2*abs(dE)/(T*rho);                                 % Time of flight
end