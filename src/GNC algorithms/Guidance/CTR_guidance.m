%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 21/05/21
% File: CTR_guidance.m 
% Issue: 0 
% Validated: 20/05/21

%% Chebyshev Trajectory Regression Guidance %%
% This script contains the function to compute the control law by means of the CTRG guidance core.

% Inputs: - scalar order, the order of the regression
%         - vector tspan, the time evolution along which the trajectory is
%           to be regressed
%         - array S, the relative state space trajectory (position and velocity) to
%           regress

% Output: - vector Cp, the position evolution trajectory
%         - vector Cv, the velocity evolution trajectory
%         - vector Cg, the acceleration evolution trajectory
%         - vector Ci, the integral of the position evolution trajectory

% New versions: 

function [Cp, Cv, Cg, Ci] = CTR_guidance(order, tspan, S)
    %Constants 
    m = 6;          %Relative phase space dimension
        
    %Divide the trajectory
    p = S(:,1:3).';                     %Position evolution
    v = S(:,4:6).';                     %Velocity evolution
    
    %Preallocation of the coefficients 
    Cp = zeros(m/2, order);             %Position regression coefficients
    Cv = zeros(m/2, order);             %Velocity regression coefficients
    
    %Regress the position 
    for i = 1:size(p,1)
        Cp(i,:) = chebyshev_coefficients(tspan, p(i,:), order);
    end
    
    %Regress the velocity
    for i = 1:size(v,1)
        Cv(i,:) = chebyshev_coefficients(tspan, v(i,:), order);
    end
    
    %Regress the acceleration 
    Cg = fcheb_derivative(Cv, tspan(end), tspan(1));

    %Regress the integral of the position
    Ci = fcheb_integral(Cp, tspan(end), tspan(1));
end