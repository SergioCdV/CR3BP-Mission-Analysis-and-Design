%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/12/24
% File: test_1.m 
% Issue: 0 
% Validated: 

%% Test 1 %%
% This scripts provides a test interface for the rest of the library
% functions

% Test 1 is concerned with the generation of the basic objects concerning a
% CR3BP system

%% Input parameters 
% System characteristics
mu = 0.0121505856;                       % Mass parameter for the Earth-Moon system

%% Basic system 
% Create a CR3BP system 
EarthMoon = src.Systems.CR3BPSystem( 'mu', mu );

%% Equilibrium points of the system 
figure 
hold on
scatter(EarthMoon.R(1,:), EarthMoon.R(2,:), 'k', 'filled');
labels = {'$M_1$', '$M_2$'};
text(EarthMoon.R(1,:), EarthMoon.R(2,:)-0.1, labels);

labels = {'$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$'};
for i = 1:size(EarthMoon.LP.r, 2)
    scatter(EarthMoon.LP.r(1,i), EarthMoon.LP.r(2,i), 'k+');
end
text(EarthMoon.LP.r(1,:), EarthMoon.LP.r(2,:)+0.1, labels);
grid on;
xlabel('$x$');
ylabel('$y$');

%% Augmented potential curves 
d = linspace(-2, +2, 100);                    % Sampling of the configuration space
[X, Y] = meshgrid(d, d);                      % Sampled configuration space

R1 = sqrt( (X + mu).^2 + Y.^2 );              % Relative position to the first primary
R2 = sqrt( (X - (1-mu)).^2 + Y.^2 );          % Relative position to the second primary

Ug = - (1 - mu) ./ R1 - mu ./ R2;             % Gravitational term
Uc = - 0.5 * (X.^2 + Y.^2);                   % Centrifugal potential
Uaug = Ug + Uc;                               % Total potential

rho = -3.5:1e-1:0;                            % Distribution of isocurves
contour(X, Y, Uaug, rho);
colorbar;

surf(X, Y, Uaug, 'FaceColor', 'red', 'EdgeColor', 'none');
camlight headlight; lighting phong
zlim([-3.5 0])
view([60 70])
zlabel('$\tilde{U}$');
alpha(0.9);

%% Zero-surface manifolds
% Compute the Zero Velocity Curve
r = src.Systems.CR3BPSystem.ZeroVelocityCurve(mu, 3.1776, false);

% Compute the Zero Velocity Surface 
r = src.Systems.CR3BPSystem.ZeroVelocitySurface(mu, 3.1776, false);

% Compute the Complementary Zero Surface 
s = src.Systems.CR3BPSystem.ComplementaryZeroSurface(mu, 3.1776, false);

%% Transformation between reference frames 
T = src.Systems.CR3BPSystem.Kepler2Synodic(mu, 1, 0, false);

% Transformation of the vector
r_earth = T * [r; ones(1,size(r,2))];

T = src.Systems.CR3BPSystem.Kepler2Synodic(mu, 1, 0, true);

% Transformation of the vector
r_syn = T * r_earth;
