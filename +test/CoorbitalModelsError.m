%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: CoorbitalModelsError.m 
% Issue: 0 
% Validated: 

%% Analysis of the different co-orbital dynamics models %%
% This scripts provides a test interface to quantify the error of the
% different models of the co-orbital CR3BP dynamics

clear;
close all;
src.graphics.set_graphics()

%% Basic system 
% Create a CR3BP system 
EarthMoon = src.Systems.EarthMoon();

%% IVP
t0 = 0;                                      % Initial clock
tf = 2*pi;                                   % Final clock

s0 = [0.8431 0 0 0 0.1874 0.4000].';         % State vector of a target vertical orbit
n = length(s0);                              % Dimensionality of the problem

TargetOrbit = src.Orbit(n, EarthMoon);
TargetOrbit.State = s0;
TargetOrbit.t = 0;

targetCR3BPIVP = src.DynamicalSystems.IVP( EarthMoon, s0, t0 );

s0 = [1.0406 0 0.1735 0 -0.0770 0].';       % State vector of a chaser butterfly orbit
s0 = [0.8431 0 0 0 0.1874 0.4001].';

ChaserOrbit = src.Orbit(n, EarthMoon);
ChaserOrbit.State = s0;
ChaserOrbit.t = 0;

chaserCR3BPIVP = src.DynamicalSystems.IVP( EarthMoon, s0, t0 );

%% Integration
% Define the solver to be used
options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14 );
integrator = src.DynamicalSystems.HybridSolver( @ode113, options );

% Integration 
dt = .001;                              % Maximum time step
tspan = [t0 tf dt];                     % Continuous horizon

Solver = integrator.configure( targetCR3BPIVP );
[t, j, tgt_y, stats] = Solver.solve( tspan );

Solver = integrator.configure( chaserCR3BPIVP );
[t, j, chs_y, stats] = Solver.solve( tspan );

if 0
    figure 
    view(3)
    hold on;
    plot3(tgt_y(1,:), tgt_y(2,:), tgt_y(3,:))
    plot3(chs_y(1,:), chs_y(2,:), chs_y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

% Relative state 
rel_state = chs_y - tgt_y;

%% Co-orbital IVP 
% Co-orbital problem
CoProblem = src.Systems.CoCR3BPSystem();

% Augment the system with the reference target
CoProblem = CoProblem .* EarthMoon;

%% Integration with the Newton model 
coCR3BPIVP = src.DynamicalSystems.IVP( CoProblem, [tgt_y(:,1); rel_state(1:6,1)], t0 );

Solver = integrator.configure( coCR3BPIVP );
[t, j, rel_y, stats] = Solver.solve( tspan );

if 0
    figure 
    view(3)
    hold on;
    plot3(tgt_y(1,:), tgt_y(2,:), tgt_y(3,:))
    plot3(rel_y(7,:) + tgt_y(1,:), rel_y(8,:) + tgt_y(2,:), rel_y(9,:) + tgt_y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

% Relative state 
dstate = chs_y - tgt_y - rel_y(7:12,:);
error(1,:) = sqrt( dot(dstate, dstate, 1) );

%% Integration with the Encke model 
coCR3BPIVP.System.params{1} = src.Systems.ModelsCoCR3BP.Encke;

Solver = integrator.configure( coCR3BPIVP );
[t, j, rel_y, stats] = Solver.solve( tspan );

if 0
    figure 
    view(3)
    hold on;
    plot3(tgt_y(1,:), tgt_y(2,:), tgt_y(3,:))
    plot3(rel_y(7,:) + tgt_y(1,:), rel_y(8,:) + tgt_y(2,:), rel_y(9,:) + tgt_y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

% Relative state 
dstate = chs_y - tgt_y - rel_y(7:12,:);
error(2,:) = sqrt( dot(dstate, dstate, 1) );

%% Integration using the first order model
coCR3BPIVP.System.params{1} = src.Systems.ModelsCoCR3BP.Linear;

Solver = integrator.configure( coCR3BPIVP );
[t, j, rel_y, stats] = Solver.solve( tspan );

if 0
    figure 
    view(3)
    hold on;
    plot3(tgt_y(1,:), tgt_y(2,:), tgt_y(3,:))
    plot3(rel_y(7,:) + tgt_y(1,:), rel_y(8,:) + tgt_y(2,:), rel_y(9,:) + tgt_y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

% Relative state 
dstate = chs_y - tgt_y - rel_y(7:12,:);
error(3,:) = sqrt( dot(dstate, dstate, 1) );

%% Integration using the second order model
coCR3BPIVP.System.params{1} = src.Systems.ModelsCoCR3BP.Order2;

Solver = integrator.configure( coCR3BPIVP );
[t, j, rel_y, stats] = Solver.solve( tspan );

if 0
    figure 
    view(3)
    hold on;
    plot3(tgt_y(1,:), tgt_y(2,:), tgt_y(3,:))
    plot3(rel_y(7,:) + tgt_y(1,:), rel_y(8,:) + tgt_y(2,:), rel_y(9,:) + tgt_y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

% Relative state 
dstate = chs_y - tgt_y - rel_y(7:12,:);
error(4,:) = sqrt( dot(dstate, dstate, 1) );

%% Integration using the third order model
coCR3BPIVP.System.params{1} = src.Systems.ModelsCoCR3BP.Order3;

Solver = integrator.configure( coCR3BPIVP );
[t, j, rel_y, stats] = Solver.solve( tspan );

if 1
    figure 
    view(3)
    hold on;
    plot3(tgt_y(1,:), tgt_y(2,:), tgt_y(3,:))
    plot3(rel_y(7,:) + tgt_y(1,:), rel_y(8,:) + tgt_y(2,:), rel_y(9,:) + tgt_y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

% Relative state 
dstate = chs_y - tgt_y - rel_y(7:12,:);
error(5,:) = sqrt( dot(dstate, dstate, 1) );

%% Integration using the Libration model
coCR3BPIVP.System.params{1} = src.Systems.ModelsCoCR3BP.Libration;

Solver = integrator.configure( coCR3BPIVP );
[t, j, rel_y, stats] = Solver.solve( tspan );

if 0
    figure 
    view(3)
    hold on;
    plot3(tgt_y(1,:), tgt_y(2,:), tgt_y(3,:))
    plot3(rel_y(7,:) + tgt_y(1,:), rel_y(8,:) + tgt_y(2,:), rel_y(9,:) + tgt_y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

% Relative state 
dstate = chs_y - tgt_y - rel_y(7:12,:);
error(6,:) = sqrt( dot(dstate, dstate, 1) );

%% Integration using Richardson's model
cn = src.Systems.CR3BPSystem.LegendreCoefficients(EarthMoon.mu, 2, EarthMoon.LP.gamma(2), 2);
coCR3BPIVP.System.params{1} = src.Systems.ModelsCoCR3BP.Richardson;
coCR3BPIVP.System.params{2} = [coCR3BPIVP.System.params{2}; cn(end)];

Solver = integrator.configure( coCR3BPIVP );
[t, j, rel_y, stats] = Solver.solve( tspan );

if 0
    figure 
    view(3)
    hold on;
    plot3(tgt_y(1,:), tgt_y(2,:), tgt_y(3,:))
    plot3(rel_y(7,:) + tgt_y(1,:), rel_y(8,:) + tgt_y(2,:), rel_y(9,:) + tgt_y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end

% Relative state 
dstate = chs_y - tgt_y - rel_y(7:12,:);
error(7,:) = sqrt( dot(dstate, dstate, 1) );

%% Integration using the N-th order model

%% Results
figure 
plot( log(t), log(error) );
legend('Newton', 'Encke', 'Linear', '2nd Order', '3rd Order', 'Libration', 'Richardson')
grid on; 
xlabel('$\textrm{log}\,t$')
ylabel('$\textrm{log}\,\mathbf{e}$')