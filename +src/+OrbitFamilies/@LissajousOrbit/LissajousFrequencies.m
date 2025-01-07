%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 06/01/25
% File: LissajousFrequencies.m 
% Issue: 0 

%% Lissajous frequencies %% 
% This functions allows to generate a Lissajous orbit seed

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - scalar L, the ID of the libration point around which the
%           Lissajous orbit is to be generated
%         - scalar gamma, the distance of the libration point to the
%           nearest primary

% Output: - vector w, the frequencies of the Lissajous orbit
%         - scalar kap, the coupling constraint in the xy plane

function [w, kap] = LissajousFrequencies(mu, L, gamma)
    % Orbit parameters (frequencies)
    cn = src.Systems.CR3BPSystem.LegendreCoefficients(mu, L, gamma, 2);     % Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(end);                                                           % Legendre coefficient c_2 (equivalent to mu)
    
    % Orbit frequencies
    w(1,1)  = sqrt( 0.5 * (2 - c2 + sqrt(9 * c2^2 - 8 * c2)) );             % In-plane frequency                 
    w(2,1)  = sqrt(c2);                                                     % Out of plane frequency
    kap = (w(1,1)^2 + 1 + 2 * c2) / (2 * w(1,1));                           % Contraint on the planar amplitude
end