%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/12/24
% File: test_2.m 
% Issue: 0 
% Validated: 

%% Test 2 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 2 is concerned with the generation of the basic objects representing
% periodic orbits in a CR3BP system

src.graphics.set_graphics()

%% Basic system 
% Create a CR3BP system 
EarthMoon = src.Systems.EarthMoon();

%% Evaluation of the Jacobian of the absolute dynamics
r = EarthMoon.ZeroVelocitySurface(EarthMoon.mu, 3.1776, false);
s = [r; zeros(3, size(r,2))];
J = EarthMoon.JacobianCR3BP(EarthMoon.mu, s);

%% Evaluation of the vector field of the absolute dynamics
s = [r; zeros(3, size(r,2))];
ds = EarthMoon.Dynamics(0, 0, s, zeros(3,size(s,2)), EarthMoon.params);

%% Solving of an IVP
s0 = [0.8431 0 0 0 0.1874 0.4000].';         % State vector of a vertical orbit
t0 = 0;                                      % Initial clock
tf = 2*pi;                                   % Final clock

%% System definition 
% Initial Value Problem 
CR3BPIVP = src.DynamicalSystems.IVP( EarthMoon, s0, t0 );

% Integrator 
options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14 );
integrator = src.DynamicalSystems.HybridSolver( @ode113, options );

%% Integration
% Configuration 
Solver = integrator.configure( CR3BPIVP );

% Integration 
dt = .001;                               % Maximum time step
tspan = [t0 tf dt];                     % Continuous horizon

[t, y, stats] = Solver.solve( tspan );

%% Results 
if 0
    figure 
    plot3(y(1,:), y(2,:), y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

%% Variational equations 
s = [y(:,1); reshape(eye(6), [], 1)];
[ds] = src.Systems.VariationalCR3BP.VariationalEquationsCR3BP(0, 0, s, [6; EarthMoon.mu]);

%% Integration of the full system 
% Create the complete variational system 
VarSystem = src.Systems.VariationalCR3BP( EarthMoon.StateDim );
CompleteSystem = EarthMoon .* VarSystem;

% Initial Value Problem 
VarCR3BPIVP = src.DynamicalSystems.IVP( CompleteSystem, [s0; reshape(eye(size(s0,1)), [], 1)], t0 );

% Configuration 
Solver = integrator.configure( VarCR3BPIVP );

% Solve the system 
[t, y, stats] = Solver.solve( tspan );

%% Results 
if 1
    figure 
    plot3(y(1,:), y(2,:), y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end