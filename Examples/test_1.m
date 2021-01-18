%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/09/20
% File: test_1.m 
% Issue: 0 
% Validated: 

%% Test 1 %%
% This scripts provides a test interface for the rest of the library
% functions. Test 1 is concerned with basic reference systems
% transformations and general analysis of a CR3BP phase space topology.

%% Test values and constants
%Set up search path 
setup_path();

%Configuration space
x = -2:1e-2:2;          %Synodic x axis
y = x;                  %Synodic y axis
z = 0;                  %Synodic plane
[X,Y] = meshgrid(x,y);  %Configuration space grid

%System characteristics
mu = 0.01215;           %Mass parameter for the Earth-Moon system
T = 24*3600;            %Synodic orbital period
d = 300e3;              %Distance between primaries
Vp = 27e3;              %Orbital velocity of the most massive primary

%Reference frame transformations 
epoch = 1;                               %Epoch at which computations are to be done
n = 2*pi/T;                              %Synodic frequency
R = [cos(n*epoch); -sin(n*epoch); 0];    %Dimensional position vector
V = [1; 0; 0];                           %Dimensional velocity vector
direction = 0;                           %Conversion direction

%Plotting parameters
rho = 1000;             %Iso-potential lines density

%% Functions
%System analysis and libration points 
L = libration_points(mu); 
C = jacobiHill_values(mu); 
c = legendre_coefficients(mu, 1, L(end,1), 10);

%Normalization procedure
epoch = dimensionalizer(d, T, Vp, epoch, 'Epoch', 1);
R = dimensionalizer(d, T, Vp, R, 'Position', 1);
V = dimensionalizer(d, T, Vp, V, 'Velocity', 1);

%Reference frames transformations  
T = inertial2synodic(epoch, [R;V], direction); 
R = inertial2synodic(epoch, R, direction);
R = synodic2lagrange(mu, L(end,1), 1, R, 0); 
R = synodic2lagrange(mu, L(end,1), 1, R, 1);

%Potential and zero-velocity curves computation
for i = 1:size(X,1)
    for j = 1:size(X,2)
        r = [X(i,j); Y(i,j); z];
        s = [r; zeros(3,1)];
        U(i,j) = augmented_potential(mu,r);
        J(i,j) = jacobi_constant(mu,s);
    end
end
figure(1)
contour(X,Y,U,rho); 
xlabel('Synodic x coordinate'); 
ylabel('Synodic y coordinate');
title('Potential function');
grid on 
colorbar

figure(2)
hold on
contour(X,Y,J,rho);
plot(L(1,1), L(2,1), 'ok');
plot(L(1,2), L(2,2), 'ok');
plot(L(1,3), L(2,3), 'ok');
plot(L(1,4), L(2,4), 'ok');
plot(L(1,5), L(2,5), 'ok');
plot(-mu, 0, 'or');
plot(1-mu, 0, 'or');
hold off
xlabel('Synodic x coordinate'); 
ylabel('Synodic y coordinate');
title('Zero velocity curves');
legend('Zero-velocity curves', 'L_1', 'L_2', 'L_3', 'L_4', 'L_5', 'M_1', 'M_2'); 
grid on 
colorbar
