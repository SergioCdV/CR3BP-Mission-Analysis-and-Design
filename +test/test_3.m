%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 25/12/24
% File: test_3.m 
% Issue: 0 
% Validated: 

%% Test 3 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 3 is concerned with the generation of the basic objects representing
% periodic orbits in a CR3BP system

src.graphics.set_graphics()

%% Basic system 
% Create a CR3BP system 
EarthMoon = src.Systems.EarthMoon();

% Create the variational system 
VarSystem = src.Systems.VariationalCR3BP( EarthMoon.StateDim );

% Create the complete variational system 
CompleteSystem = EarthMoon .* VarSystem;

%% Solving of an IVP
s0 = [0.8431 0 0 0 0.1874 0.4000].';         % State vector of a vertical orbit
n = length(s0);                              % Dimensionality of the problem

myOrbit = src.Systems.Orbit(n, EarthMoon);

STM = src.Systems.STM( n );                  % Initial conditions of the STM
Phi = reshape(STM.Phi, [], 1);              

t0 = 0;                                      % Initial clock
tf = 2*pi;                                   % Final clock

% Initial Value Problem 
VarCR3BPIVP = src.DynamicalSystems.IVP( CompleteSystem, [s0; Phi], t0 );

% Integrator 
options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14 );
integrator = src.DynamicalSystems.HybridSolver( @ode113, options );

% Configuration 
Solver = integrator.configure( VarCR3BPIVP );

% Integration 
dt = .001;                              % Maximum time step
tspan = [t0 tf dt];                     % Continuous horizon

%% Integration
% Solve the system 
[t, y, stats] = Solver.solve( tspan );

if 1
    figure 
    plot3(y(1,:), y(2,:), y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

%% Tests on the methods of the Orbit object 
% Add the phase space trajectory to the Orbit 
trajectory{1} = t; 
trajectory{2} = y(1:n,:);
myOrbit.State = trajectory;

% Add the STM to the Orbit 
myOrbit.STM = STM;

%% Tests on the methods of the STM object 
STM.Phi = y(7:end,:);
CGT = STM.CauchyGreenTensor( STM.Phi );

% Symplecticity error 
error = src.CheckSymplecticity( STM.Phi(:,7:12) );

if 0
    figure 
    plot(t, STM.StabilityIndex(1,:)); 
    grid on;
    xlabel('$t$')
    ylabel('$s$')
end