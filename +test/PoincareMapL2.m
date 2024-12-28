%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/12/24
% File: PoincareMapL2.m 
% Issue: 0 
% Validated: 

%% Poincaré Map in the vecinity of L1 %%
% This scripts provides a test interface for the PoincareMap class by
% computing the 1st-return map to the section y = 0 for the vecinity of the
% L1 point in the Earth Moon system 

clear;
close all;
src.graphics.set_graphics()

%% Basic system 
% Create a CR3BP system 
EarthMoon = src.Systems.EarthMoon();

%% Meshing the vecinity of L2
Hlevel = EarthMoon.LP.J(1);                  % Energy level of the map

% Blob of phase space
s0 = [0.8431 0 0 0 0.1874 0.4000].';         % State vector of a vertical orbit
n = length(s0);                              % Dimensionality of the problem

NumOrbits = 100;                             % Number of initial conditions to analyze
Blob = {};                                   % Pre-allocation
l = 1; 

x = EarthMoon.LP.r(1,1) + linspace(-0.001, 0.001, NumOrbits);
y = EarthMoon.LP.r(2,1) + linspace(-0.05, 0.05, NumOrbits);
z = EarthMoon.LP.r(3,1) + linspace(-0.05, 0.05, NumOrbits);

[X, Y, Z] = meshgrid(x, y, z);

for i = 1:NumOrbits
    for j = 1:NumOrbits
        for k = 1:NumOrbits
            % Sampling 
            x = X(i,j,k); 
            y = Y(i,j,k);
            z = Z(i,j,k);
                    
            if ( randi([0 1]) == 0 )
                signo = 1; 
            else
                signo = -1;
            end

            vx = signo * sqrt( 2 * (Hlevel - EarthMoon.JacobiConstant(EarthMoon.mu, [x y z 0 0 0].')) );

            if ( imag(vx) ~= 0 )
                break;
            else
                s0 = [x y z vx 0 0].';
            
                % Definition of the objects
                myBlobOrbit = src.Orbit(n, EarthMoon);
                myBlobOrbit.State = s0;
                myBlobOrbit.t = 0;
                Blob{l} = myBlobOrbit;
                l = l + 1;
            end
        end
    end
end

%% Create the Poincare map object 
mySoS = @(t, j, s, u, params)( s(2,:) );       % Transversal surface of section of the flow
myNumberReturns = 1;                           % Number of returns allowed by the map

myMap = src.PoincareMap(myNumberReturns, mySoS);

%% Compute the map
% Define the solver to be used
options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14 );
Solver = src.DynamicalSystems.HybridSolver( @ode113, options );

% Compute the return map for the selected phase space blob
[r, OrbitCollection] = myMap.Compute( Solver, Blob );                                      

%% Plot results 
if 1 && ~isempty(r)
    figure
    scatter(r(1,:), r(3,:), '.', 'k')
    xlabel('$x$');
    ylabel('$z$');
    grid on;
end