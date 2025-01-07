%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/24
% File: OrbitSeed.m 
% Issue: 0 

%% Lissajous Orbit seed %% 
% This functions allows to generate a Lissajous orbit seed

% Inputs: - obj Orbit, the Lissajous orbit of interest
%         - vector Amp, the amplitudes of the Lissajous orbit
%         - vector theta[2 x N], the set of phases of the orbit 
%         - vector freq, the set of frequencies of the orbit 
%         - double kap, the amplitude constraint in the xy plane

% Output: - array seed [6 x N], containing the required initial solution seed

function [seed] = OrbitSeed(obj, Amp, theta, freq, kap) 
    % Sanity checks 
    if ( ~exist("freq", "var") )
        freq = obj.OrbitFrequencies;
    end

    if ( ~exist("kap", "var") )
        kap = obj.kap;
    end

    if ( ~exist("theta", "var") )
        theta = zeros(2,1);
    end

    if ( length(Amp) < 2 )
        O = zeros(1, 2-size(Amp,1) );
        Amp = [Amp O];
    end

    % Parameters of the Lyapunov orbit
    Ax = Amp(1);                            % In-plane trajectory
    Az = Amp(2);                            % Out-of-plane trajectory

    phi = theta(1,:);                       % In-plane phase
    psi = theta(2,:);                       % Out-of-plane phase

    % Seed trajectory
    x = -1.0 * Ax * cos(phi);               % X relative coordinate
    y =  kap * Ax * sin(phi);               % Y relative coordinate
    z =        Az * sin(psi);               % Z relative coordinate
    vx = freq(1)  * Ax * sin(phi);          % Vx relative velocity
    vy = kap * freq(1) * Ax * cos(phi);     % Vy relative velocity
    vz = freq(2) * Az * cos(psi);           % Vz relative velocity  
        
    % Output seed
    seed = [x; y; z; vx; vy; vz];
end