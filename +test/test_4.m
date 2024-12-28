%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 28/12/24
% File: test_4.m 
% Issue: 0 
% Validated: 

%% Test 4 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 4 is concerned with the solving of a co-orbital IVP with variational
% augmentation

src.graphics.set_graphics()

%% Basic, reference system 
% Create a CR3BP system 
EarthMoon = src.Systems.EarthMoon();

% Create the variational system 
VarSystem = src.Systems.VariationalCR3BP( EarthMoon.StateDim );

% Create the complete variational system 
TargetCompleteSystem = EarthMoon .* VarSystem;

%% Co-orbital system 
% Create a CR3BP system 
CoOrbital = src.Systems.CoCR3BPSystem();

% Create the variational system 
CoVarSystem = src.Systems.VariationalCoCR3BP( 6 );

% Create the complete variational system 
CoOrbitalCompleteSystem = CoOrbital .* CoVarSystem;

%% Complete system 
CoOrbitalCompleteSystem = CoOrbitalCompleteSystem .* TargetCompleteSystem;

%% Solving of the IVP
s0 = [0.8431 0 0 0 0.1874 0.4000].';         % State vector of a vertical orbit
n = length(s0);                              % Dimensionality of the problem

myOrbit = src.Orbit(n, EarthMoon);           % Basic orbit

STM = src.STM( n );                          % Initial conditions of the STM
Phi = reshape(STM.Phi, [], 1);              

t0 = 0;                                      % Initial clock
tf = 2*pi;                                   % Final clock

dstate = zeros(6,1);                         % Relative state

% Initial Value Problem 
CoCR3BPIVP = src.DynamicalSystems.IVP( CoOrbitalCompleteSystem, [s0; Phi; dstate; Phi], t0 );

% Integrator 
options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14 );
integrator = src.DynamicalSystems.HybridSolver( @ode113, options );

% Configuration 
Solver = integrator.configure( CoCR3BPIVP );

% Integration 
dt = .001;                              % Maximum time step
tspan = [t0 tf dt];                     % Continuous horizon

%% Integration
% Solve the system 
[t, j, y, stats] = Solver.solve( tspan );

%% Results
% Chaser motion 
tgt_y = y(1:EarthMoon.StateDim,:);
chs_y = tgt_y + y(CoOrbitalCompleteSystem.OriginalStateDim(1) + 1:CoOrbitalCompleteSystem.OriginalStateDim(1) + CoOrbitalCompleteSystem.PhaseSpaceDim(2),:);

if 1
    figure 
    view(3)
    hold on
    plot3(tgt_y(1,:), tgt_y(2,:), tgt_y(3,:))
    plot3(chs_y(1,:), chs_y(2,:), chs_y(3,:))
    grid on; 
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
end