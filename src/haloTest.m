%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/09/20
% File: haloTest.m 
% Issue: 0 
% Validated: 

%% Test Halo %%
% This scripts provides a useful interface to test periodic orbits computations and
% their associated manifolds

%% General setup 
clc
format long;

setup_path(); 

%% Initial conditions 
%Gravitational parameter 
mu =  3.054248395728148e-06;

%Libration points calculation 
L = libration_points(mu); 

%Selection of the libration point
libration_point = 'L2'; 
switch (libration_point)
    case 'L1'
        Lx = L(1,1);
    case 'L2' 
        Lx = L(1,2);
    case 'L3' 
        Lx = L(1,3);
    case 'L4' 
        Lx = L(1,4);
    case 'L5' 
        Lx = L(1,5);
    otherwise
        disp('No valid options was selected.'); 
        return;
end

%Amplitude definition 
Ax = 1e-4;  %XY plane amplitude
Az = 0;     %Z amplitude

%Initial conditions computation (first order analytical approximation)
bmu = mu*abs(Lx-1+mu)^(-3)+(1-mu)*abs(Lx+mu)^(-3);  %Reduced gravitational parameter
nu = sqrt(-0.5*(bmu-2-sqrt(9*bmu^2-8*bmu)));
tau = -(nu^2+2*bmu+1)/(2*nu);
x0 = Lx-Ax;                                         %Initial position 
vy0 = -Ax*tau*nu;                                   %Initial velocity
phi0 = eye(6);                                      %Initial STM
phi0 = reshape(phi0, [1 36]);                       %Initial STM
x0 = [x0 0 0 0 vy0 0 phi0];                         %Initial conditions

%% Differential correction targetting 
vTol = 1e-11;           %Velocity tolerance
gTol = 1e-12;           %General tolerance
iterMax = 10;           %Maximum number of iterations
go_on = true;           %Convergence flag
iter = 1;               %Initial iteration

while ((go_on) && (iter < iterMax))
    %Dynamics integration 
    [sol, phi, elaps_time] = integration_loop(mu, x0);
    
    %Targetting and convergence
    if (abs(sol(end,4)) < vTol)
        go_on = false;
    else
        %Differential correction
        dVy = sol(end,4)/(phi(4,5)-(sol(end,4)/sol(end,5))*phi(2,5));   %Required change in velocity
        x0(5) = x0(5)-dVy;                                              %New initial conditions
        
        %Display iteration 
        fprintf('Differential correction scheme iteration: %i \n', iter);
        
        %Update iteration
        iter = iter+1;
    end
end

%Integration of the whole periodic orbit 
RelTol = 2.25e-14; 
AbsTol = 1e-20;
options = odeset('RelTol', RelTol, 'AbsTol', AbsTol);

n = max(roots([1 1 -size(x0,2)]));  %Degrees of freedom
dt = 1e-4;                          %Time step
tspan = 0:dt:2*elaps_time;          %Timespan

dir = 1;                            %Temporal direction of the flow
flagVar = true;                     %Restriction flag for the number of degrees of freedom
if (~flagVar)
    x0 = x0(end,1:n);
end

%Integration
[~, r] = ode113(@(t,x)cr3bp_equations(mu, dir, flagVar, t, x), tspan, x0, options);

%Plotting the orbit
plot3(r(:,1), r(:,2), r(:,3))

%% Invariant manifolds
%Manifold computation
rho = 50;
manifold = ['U' 'U' 'S' 'S'];
branch = ['L' 'R' 'L' 'R'];
leftUM =  invariant_manifold(mu, manifold(1), branch(2), r, rho, 0:1e-3:5*elaps_time);
%rightUM = invariant_manifold(mu, manifold(2), branch(1), r, rho, elaps_time);
%leftSM =  invariant_manifold(mu, manifold(3), branch(1), r, rho, elaps_time);
%rightSM = invariant_manifold(mu, manifold(4), branch(1), r, rho, 0elaps_time);

%Plotting and results 
figure(2)
hold on 
for i = 1:size(leftUM,1)
    auxM1 = shiftdim(leftUM(i,:,:));
%     auxM2 = shiftdim(rightUM(i,:,:));
%     auxM3 = shiftdim(leftSM(i,:,:));
%     auxM4 = shiftdim(rightSM(i,:,:));
    
    plot3(auxM1(:,1), auxM1(:,2), auxM1(:,3), 'b');
%     plot3(auxM2(:,1), auxM2(:,2), auxM2(:,3));
%     plot3(auxM3(:,1), auxM3(:,2), auxM3(:,3));
%     plot3(auxM4(:,1), auxM4(:,2), auxM4(:,3));
end
plot3(r(:,1), r(:,2), r(:,3), 'r');
hold off
grid on
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');

%% Family of orbits
%Pseudo-arclength continuation: to be implemented.

%% Auxiliary functions 
function [Sol, Phi, elaps] = integration_loop(mu, x0)
    %Integration options
    RelTol = 2.5e-14;
    AbsTol = 1e-22;
    options = odeset('RelTol', RelTol, 'AbsTol', AbsTol, 'Events', @(t,x)x_crossing(t,x)); 
    
    %Integration parameters 
    n = 6;                  %Phase space dimension
    dt = 1e-3;              %Time step
    dir = 1;                %Temporal direction of the flow
    flagVar = true;         %Restriction flag of the number of degrees of freedom to integrate
    
    %Main integration loop
    %Integration span set upt
    Dt = 0:dt:1e3;
        
    %Dynamics integration 
    [~,x_aux,~,~,~] = ode113(@(t,x)cr3bp_equations(mu, dir, flagVar, t, x), ...
                             Dt, x0(end,:), options);
    
    %Output definition of the STM
    elaps = dt*size(x_aux,1);
    Sol = x_aux;
    Phi = reshape(Sol(end,n+1:end), [n n]);
end
