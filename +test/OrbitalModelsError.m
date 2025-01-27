%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 05/01/25
% File: OrbitalModelsError.m 
% Issue: 0 
% Validated: 

%% Analysis of the different orbital dynamics models %%
% This scripts provides a test interface to quantify the error of the
% different models of the CR3BP dynamics

clear;
close all;
src.graphics.set_graphics()

%% Basic system 
% Create a CR3BP system 
EarthMoon = src.Systems.EarthMoon();

%% IVP
t0 = 0;                                      % Initial clock
tf = pi;                                   % Final clock

s0 = [EarthMoon.LP.r(1,1) 0 0 0 1E-7 0].';         % State vector of a target vertical orbit

s0 = [0.824024728136525 0 -0.054501847320725 0 0.164671964079122 0].';
n = length(s0);                              % Dimensionality of the problem

TargetOrbit = src.Orbit(n, EarthMoon);
TargetOrbit.State = s0;
TargetOrbit.t = 0;

CR3BPIVP = src.DynamicalSystems.IVP( EarthMoon, s0, t0 );

options = odeset('AbsTol', 1E-10, 'RelTol', 2.25E-10 );
integrator = src.DynamicalSystems.HybridSolver( @ode113, options );

% Integration 
dt = .001;                              % Maximum time step
tspan = [t0 tf dt];                     % Continuous horizon

%% Integration with Newton model
CR3BPIVP.System.params{2} = [CR3BPIVP.System.params{2}; EarthMoon.LP.gamma(1)];

Solver = integrator.configure( CR3BPIVP );
[t, j, Newton_y, stats] = Solver.solve( tspan );

%% Integration with the Encke model 
CR3BPIVP.System.params{1} = src.Systems.ModelsCR3BP.Encke;

Solver = integrator.configure( CR3BPIVP );
[t, j, Encke_y, stats] = Solver.solve( tspan );

% Relative state 
dstate = Newton_y - Encke_y;
error(1,:) = sqrt( dot(dstate, dstate, 1) );

%% Integration using the N-th order model, as a co-orbital propagation from a given libration point
% Model parameters
Lidx = 1;                       
CR3BPIVP.System.params{1} = src.Systems.ModelsCR3BP.OrderN;
CR3BPIVP.System.params{2} = [EarthMoon.mu; Lidx; EarthMoon.LP.gamma(Lidx); 100];

% Co-orbital problem
CR3BPIVP.IC(1:3,:) = CR3BPIVP.IC(1:3,:) - EarthMoon.LP.r(1:3,Lidx);
CR3BPIVP.IC(1:6,:) = CR3BPIVP.IC(1:6,:) / EarthMoon.LP.gamma(Lidx);

% Integration
Solver = integrator.configure( CR3BPIVP );
[t, j, Rich_y, stats] = Solver.solve( tspan );
 
% Re-scaling and traslation of the origin
Rich_y = Rich_y * EarthMoon.LP.gamma(Lidx);
Rich_y(1:3,:) = Rich_y(1:3,:) + EarthMoon.LP.r(1:3,Lidx);

% Relative state 
dstate = Newton_y - Rich_y;
error(2,:) = sqrt( dot(dstate, dstate, 1) );

%% Results
if 1
    figure 
    view(3)
    hold on;
    plot3(Newton_y(1,:), Newton_y(2,:), Newton_y(3,:), 'b')
    plot3(Encke_y(1,:), Encke_y(2,:), Encke_y(3,:), 'r')
    plot3(Rich_y(1,:), Rich_y(2,:), Rich_y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

figure 
plot( log(t), log( error ) );
legend('Encke', 'N-th order')
grid on; 
xlabel('$\textrm{log}\,t$')
ylabel('$\textrm{log}\,\mathbf{e}$')