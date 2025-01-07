%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/01/25
% File: TimeLaw.m 
% Issue: 0 

%% Time law %% 
% This functions allows to compute the phase associated to a given epoch on
% the orbit

% Inputs: - obj Orbit, the Lissajous orbit of interest
%         - vector tspan [1 x N], the epochs at which to evaluate the orbit
%         - vector theta, the set of initial phases on the orbit 

% Output: - array phi [2 x N], containing the phases on the orbit
%         - num_rev [1 x N], containing the number of completed revolutions
%           of the orbit

function [phi, num_rev] = TimeLaw(obj, tspan, theta) 
    % Sanity checks 
    if ( ~exist("theta", "var") )
        theta = zeros(2,1);
    end

    freq = obj.OrbitFrequencies;

    % Phase on the orbit
    phi = freq .* tspan + theta;

    % Number of revolutions 
    T = 2*pi / freq(1);                     % Planar period
    num_rev = floor( tspan / T );           % Number of revolutions
end