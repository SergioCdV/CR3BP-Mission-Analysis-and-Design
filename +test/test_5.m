%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 07/01/25
% File: test_5.m 
% Issue: 0 
% Validated: 

%% Test 5 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 5 is concerned with the generation of different orbit families

src.graphics.set_graphics()

%% Basic, reference system 
% Create a CR3BP system 
EarthMoon = src.Systems.EarthMoon();

%% Lissajous orbit
% Create the orbit
LissajousOrbit = src.OrbitFamilies.LissajousOrbit( EarthMoon, 2 );

% Define the amplitudes
LissajousOrbit.OrbitAmplitudes(1:2) = [0.01, 0.02];        % Normalized units

% Generate the orbit seed 
theta = linspace(0, 16*pi, 1000);
theta = LissajousOrbit.TimeLaw( theta );
seed = LissajousOrbit.OrbitSeed( LissajousOrbit.OrbitAmplitudes, theta );
seed(1:3,:) = seed(1:3,:) + LissajousOrbit.Origin;

if 0
    figure 
    view(3)
    hold on
    labels = {'$L_2$'};
    scatter(EarthMoon.LP.r(1,2), EarthMoon.LP.r(2,2), 'k+');
    text(EarthMoon.LP.r(1,2), EarthMoon.LP.r(2,2)+0.1, labels);
    plot3(seed(1,:), seed(2,:), seed(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

%% Lyapunov orbit
% Create the orbit
LyapunovOrbit = src.OrbitFamilies.LyapunovOrbit( EarthMoon, 2 );

% Define the amplitudes
LyapunovOrbit.OrbitAmplitudes = [0.01];        % Normalized units

% Generate the orbit seed 
theta = linspace(0, 16*pi, 1000);
theta = LyapunovOrbit.TimeLaw( theta );
seed = LyapunovOrbit.OrbitSeed( LyapunovOrbit.OrbitAmplitudes, theta );
seed(1:3,:) = seed(1:3,:) + LyapunovOrbit.Origin;

if 0
    figure 
    view(3)
    hold on
    labels = {'$L_2$'};
    scatter(EarthMoon.LP.r(1,2), EarthMoon.LP.r(2,2), 'k+');
    text(EarthMoon.LP.r(1,2), EarthMoon.LP.r(2,2)+0.1, labels);
    plot3(seed(1,:), seed(2,:), seed(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

%% Halo orbit
% Create the orbit
HaloOrbit = src.OrbitFamilies.HaloOrbit( EarthMoon, 2, src.OrbitFamilies.HaloBranches.Northern );

% Define the amplitudes
HaloOrbit.OrbitAmplitudes = [0.01, 0.02];        % Normalized units

% Generate the orbit seed 
theta = linspace(0, 16*pi, 1000);
theta = HaloOrbit.TimeLaw( theta );
seed = HaloOrbit.OrbitSeed( HaloOrbit.OrbitAmplitudes, theta );
seed(1:3,:) = seed(1:3,:) + HaloOrbit.Origin;

if 1
    figure 
    view(3)
    hold on
    labels = {'$L_2$'};
    scatter(EarthMoon.LP.r(1,2), EarthMoon.LP.r(2,2), 'k+');
    text(EarthMoon.LP.r(1,2), EarthMoon.LP.r(2,2)+0.1, labels);
    plot3(seed(1,:), seed(2,:), seed(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end