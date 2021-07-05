%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/02/21
% File: polynomials_graphics.m 
% Issue: 0 
% Validated: 

%% Polynomials graphics %% 
% This script provides a function to plot several Legendre and Chebyshev
% polynomials 

%% Main computation %%
set_graphics(); 

%Initial data
theta = 0:1e-2:2*pi;             %Normalized domain
x0 = [1; 0; 0];                  %Center of the orbit

%Coorbital motion solution 
Se = [cos(theta); sin(theta); zeros(1,length(theta))] + repmat(x0, 1, length(theta));

%Libration solution
R = 1e-1; 
r = [cos(theta(500)) -sin(theta(500)) 0; sin(theta(500)) cos(theta(500)) 0; 0 0 1];
Sl = R*r*[sin(theta); zeros(1,length(theta)); cos(theta)] + repmat(Se(:,500), 1, length(theta));

%% Plots and results %% 
figure(1) 
view(3)
hold on 
plot3(Se(1,:), Se(2,:), Se(3,:));
plot3(Sl(1,:), Sl(2,:), Sl(3,:));
hold off
grid on; 
title('Continuum equilibria in two-body relative motion')
xlabel('Normalized inertial $x/a_t$ coordinate');
ylabel('Normalized inertial $y/a_t$ coordinate');
zlabel('Normalized inertial $z/a_t$ coordinate'); 
axis([0 2.5 -1.1 1.1 -0.5 0.5])
legend('Co-orbital motion', 'Libration motion', 'Location', 'northeast');

