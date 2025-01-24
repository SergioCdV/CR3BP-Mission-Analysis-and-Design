%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/01/25
% File: test_6.m 
% Issue: 0 
% Validated: 

%% Test 6 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 6 is concerned with the generation of different orbits through
% differential correction

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
seed = LissajousOrbit.OrbitSeed( LissajousOrbit.OrbitAmplitudes );
seed(1:3,:) = seed(1:3,:) + LissajousOrbit.Origin;
LissajousOrbit.State = seed;
LissajousOrbit.t = 0;

%% Differential corrector 
% Create the corrector
myCorrector = src.Correctors.PlanarCorrector();

Config.RelTol = 1E-5; 
Config.AbsTol = 1E-10;
Config.MaxIter = 100; 

myCorrector = myCorrector.Configure( Config );

% Differential correction
[LissajousOrbit, Stats] = myCorrector.SingleShootSolve( LissajousOrbit );

if 1
    figure 
    plot3(LissajousOrbit.State(1,:), LissajousOrbit.State(2,:), LissajousOrbit.State(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

%% Halo orbit
% Create the orbit
HaloOrbit = src.OrbitFamilies.HaloOrbit( EarthMoon, 2, src.OrbitFamilies.HaloBranches.Northern );

% Define the amplitudes
HaloOrbit.OrbitAmplitudes = [0.1, 0.02];        % Normalized units

% Generate the orbit seed 
seed = HaloOrbit.OrbitSeed( HaloOrbit.OrbitAmplitudes );
seed(1:3,:) = seed(1:3,:) + HaloOrbit.Origin;
HaloOrbit.State = seed;
HaloOrbit.t = 0;

% Create the corrector
myCorrector = src.Correctors.PlaneCorrector();

Config.RelTol = 1E-5; 
Config.AbsTol = 1E-10;
Config.MaxIter = 100; 

myCorrector = myCorrector.Configure( Config );

% Differential correction
[HaloOrbit, Stats] = myCorrector.SingleShootSolve( HaloOrbit );

if 1
    figure 
    plot3(HaloOrbit.State(1,:), HaloOrbit.State(2,:), HaloOrbit.State(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end
