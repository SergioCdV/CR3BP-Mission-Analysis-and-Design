%% Differential correction example: free fall motion %% 
% Sergio Cuevas del Valle % 
% December 2020 % 

%% Description %%
% This file provides a differential correction scheme applied to
% projectile motion. Initial conditions, chosen randomly, are modified to
% accomplish surpassing a certain height at a certain range. A second
% scenary aims to target orbital motion from random initial conditions.

%% Scenario 1: surpassing height at known range, varying initial velocities 
% Results: time can be fixed by dynamical constraints, and should not be corrected by brute force. Modifying the 
% flight path angle and the velocity modulus seem unstable.

%Constants 
g = 1;                          %Gravity acceleration, normalized units.

% Initial conditions 
[dt, x0, xf] = initial_conditions(); 

% Generate initial trajectory for plotting 
tf = xf(1)/x0(1);           	%End time, when the obstacle is reached 
t = 0:dt:tf;                    %Integration time
vx = x0(1)*ones(1,length(t));   %Horizontal velocity
vy = x0(2)-g*t;                 %Vertical velocity
x = x0(1)*t;                    %X trayectory
y = vy.*t-(1/2)*g*t.^2;         %Y trajectory
s1 = [x; y; vx; vy];            %State trayectory
s = s1;

% Differential correction scheme 
tol = 1e-10;        %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 1e3;      %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(2, maxIter);   

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = s(1:2,end)-xf; 
    
    %Evaluate the STM at the final state 
    STM = diag([tf tf]);
    
    %Compute the needed correction 
    dx0(:,iter) = STM\e;
   
    %Generate new trajectory 
    x0 = x0-dx0(1:2,iter);                      %New initial conditions
    tf = xf(1)/x0(1);                           %New end time, when the obstacle is reached 
    t = 0:dt:tf;                                %Integration time
    vx = x0(1)*ones(1,length(t));               %Horizontal velocity
    vy = x0(2)-g*t;                             %Vertical velocity
    x = x0(1)*t;                                %X trayectory
    y = vy.*t-(1/2)*g*t.^2;                     %Y trajectory
    s = [x; y; vx; vy];                         %State trayectory
    
    %Convergence 
    if (norm(e) < tol)
        GoOn = false;
    else
        iter = iter+1;
    end
end

% Results
if (~GoOn)
    disp('Case study 1: PASS.');
else
    disp('Case study 2: FAIL. Error convergence was not met.');
end

% Plotting
figure(1) 
hold on
plot(s1(1,:), s1(2,:), '-.r');
plot(s(1,:), s(2,:), 'b');
plot(xf(1), xf(2), 'o');
hold off 
grid on;
xlabel('Range'); 
ylabel('Altitude');
title('Differential correction results for free fall motion'); 
legend('Initial trayectory', 'Final trayectory', 'Target');

%% Scenario 2: targeting orbital motion varying both location and speed (location should not matter)
% Results: reparametrization needed. The parabolic flight equations do
% not generate Keplerian motion. Curvature has to be taken into account.
% Selecting the final conditions randomly result in poor convergence.
% The obtained curve is though parabolic. Another main problem is that the
% end time is unknown, and should be input as another
% to-be-corrected-variable. It could be interesting to compute the STM
% using the analytical expression for hyperbolic orbits.

%Constants 
mu = 1;     %Gravity acceleration, normalized units.

% Initial conditions 
[dt, x0] = initial_conditions2();   %Initial conditions
xf = [2*pi; 0.5; 0; 0];             %Orbital conditions
Phi = eye(4);                       %Initial state transition matrix
Phi = reshape(Phi, [16,1]);
x0 = [x0; Phi];                     %Augmented initial conditions

% Generate initial trajectory
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);           %Integration tolerances
tf = 1;                                                     %Desired orbital period
tspan = 0:dt:tf;                                            %Integration time
[~, s] = ode45(@(t,s)dynamics(t,s), tspan, x0, options);    %Trajectory integration

% Differential correction scheme 
tol = 1e-10;        %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 1e3;      %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(1, maxIter);   

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = s(end,2:4).'-xf(2:4); 
    
    %Evaluate the STM at the final state 
    STM = reshape(s(end, 5:end), [4,4]);
    STM = STM(2:end,3);
    
    %Compute the needed correction 
    dx0(iter) = pinv(STM)*e;
   
    %Generate new trajectory 
    x0(3) = x0(3)-dx0(iter);                                    %New initial conditions
    [~, s] = ode45(@(t,s)dynamics(t,s), tspan, x0, options);    %Trajectory integration
    
    %Convergence 
    if (norm(e) < tol)
        GoOn = false;
    else
        iter = iter+1;
    end
end

% Results
if (~GoOn)
    disp('Case study 2: PASS.');
else
    disp('Case study 2: FAIL. Error convergence was not met.');
end

% Plotting 
figure(2) 
plot(s(:,1), s(:,2));
xlabel('Range');
ylabel('Altitude');
grid on;
title('Last correction to generate a Keplerian parabolic trajectory')

%% Auxiliary functions
%Generate initial conditions
function [dt, x0, xf] = initial_conditions()
    dt = 1e-5;                  %Time step
    x0 = rand(2,1);             %Initial conditions randomly chosen in normalized coordinates for velocity
    xf = x0+abs(rand(2,1));     %Obstacle situation and heigth
end

function [dt, x0] = initial_conditions2()
    dt = 1e-3;                      %Time step
    x0 = sqrt(2/0.6);               %Initial velocity 
    x0 = [0; 0; x0; deg2rad(89)];   %Initial conditions
end

%Dynamics of the orbital problem
function [dX] = dynamics(t, s)
    %Variables de estado 
    x = s(1);       %Horizontal coordinate
    y = s(2);       %Vertical coordinate
    V = s(3);       %Velocity modulus
    gamma = s(4);   %Flight path angle
    
    Phi0 = reshape(s(5:end), [4,4]);
    
    %Gravity field
    g = (1/(1+y)^2);
    
    %Dynamics vector field
    dS = [(1/(1+y))*V*cos(gamma);
                    V*sin(gamma); 
                   -g*sin(gamma); 
       -(g/V-V/(1+y))*cos(gamma)];
   
   %Variational equations
   J = [0 -(V/(1+y)^2)*cos(gamma) (1/(1+y))*cos(gamma) -(V/(1+y))*sin(gamma);
        0                        0           sin(gamma)           V*cos(gamma); 
        0                   -2*g/(1+y)               0           -g*cos(gamma); 
        0  -(-2*g/(V*(1+y))+V/(1+y)^2)*cos(gamma) -(-g/V^2-1/(1+y))*cos(gamma) -(g/V-V/(1+y))*sin(gamma)];
    
   dPhi = J*Phi0; 
   dPhi = reshape(dPhi, [16,1]);
   
   %Complete vector field 
   dX = [dS; dPhi];
end