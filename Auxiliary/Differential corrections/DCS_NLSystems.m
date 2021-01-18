%% Differential correction example: periodic orbit in non linear systems %% 
% Sergio Cuevas del Valle % 
% December 2020 % 

%% Description %%
% This file provides a differential correction scheme to generate 
% periodic orbits in several non linear systems. Initial conditions, 
% chosen randomly, are modified to generate periodic orbit of a desired 
% period. For more details on the topic, reference may be found in Wiesel and
% Strogatz. 

%% Scenario 1: periodic orbits for the Duffing oscillator
% Results: really slow convergence for the complete model... Integration of
% the STM really adds computational effort. Convergence is achieved though for the 
% most cases.

% Set up integration
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);

% Augmented initial conditions 
[dt, x0, param, T] = initial_conditions(); 
Phi = eye(2);                       %Initial STM
Phi = reshape(Phi, [4 1]);
x0 = [x0; Phi];                     %Initial conditions

% Generate initial trajectory 
tspan = 0:dt:T; 
[~, x1] = ode113(@(t,x)duffing_dynamics(t, x, param), tspan, x0, options);
s = x1;

% Differential correction scheme 
tol = 1e-5;         %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 10;       %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(2, maxIter);

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = (s(end,1:2)-s(1,1:2)).'; 
    
    %Evaluate the STM at the final state 
    STM = reshape(s(end,3:end), 2, 2); 
    
    %Compute the needed correction 
    dx0(:,iter) = (STM-eye(2))\e;
    
    %Generate new trajectory 
    x0(1:2) = x0(1:2)-dx0(1:2,iter);                                           %New initial conditions
    [~, s] = ode113(@(t,x)duffing_dynamics(t, x, param), tspan, x0, options);  %Dynamics integration
    
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
    disp('Case study 1: FAIL. Error convergence was not met.');
end
 
% Plotting
figure(1) 
hold on
plot(x1(:,1), x1(:,2), '-.r');
plot(s(:,1), s(:,2), 'b');
hold off 
xlabel('Position'); 
ylabel('Velocity');
title('Differential correction results in the phase space'); 
legend('Initial trayectory', 'Final trayectory');

%% Scenario 2: periodic orbit in the Lorenz system 
% Results: the ordinary Newton method is singular. Initial conditions must
% be refined, as sligth perturbations are quickly constrained to the
% attractor.

% Set up integration
options = odeset('RelTol', 2.25e-10, 'AbsTol', 1e-22);

% Augmented initial conditions 
[dt, x0, param, T] = initial_conditions3(); 
Phi = eye(3);                               %Initial STM
Phi = reshape(Phi, [9 1]);
x0 = [x0; Phi];                             %Initial conditions

% Generate initial trajectory 
tspan = 0:dt:T; 
[~, x1] = ode113(@(t,x)lorenz_dynamics(t, x, param), tspan, x0, options);
x0 = [x1(end,1:3).'; Phi];
[~, x1] = ode113(@(t,x)lorenz_dynamics(t, x, param), tspan, x0, options);
s = x1;

% Differential correction scheme 
tol = 1e-5;         %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 200;       %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(4, maxIter);

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = (s(end,1:3)-s(1,1:3)).'; 
    
    %Evaluate the STM at the final state
    F1 = lorenz_dynamics(0, s(end,:), param);
    F2 = lorenz_dynamics(0, s(1,:), param);
    STM = reshape(s(end,4:end), 3, 3);
    A = [STM-eye(3) F1(1:3); F2(1:3).' 0];
    
    %Compute the needed correction 
    dx0(:,iter) = A.'*(A*A.')^(-1)*[e; 0];
    
    %Generate new trajectory 
    x0(1:3) = x0(1:3)-dx0(1:3,iter);                                          %New initial conditions  
    T = T-dx0(end,iter);                                                      %Update the period
    tspan = 0:dt:T;
    [~, s] = ode113(@(t,x)lorenz_dynamics(t, x, param), tspan, x0, options);  %Dynamics integration
    
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
hold on 
plot3(s(end,1), s(end,2), s(end,3), 'ob');
plot3(s(:,1), s(:,2), s(:,3), 'r');
plot3(s(1,1), s(1,2), s(1,3), 'ok');
hold off
xlabel('x coordinate'); 
ylabel('y coordinate');
zlabel('z coordinate')
title('Differential correction results in the phase space'); 
legend('Final state', 'Trayectory', 'Initial conditions');
grid on

%% Scenario 3: periodic orbits for the Lotka-Volterra model
% Results: convergence.

% Set up integration
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);

% Augmented initial conditions 
[dt, x0, param, T] = initial_conditions4(); 
Phi = eye(2);                       %Initial STM
Phi = reshape(Phi, [4 1]);
x0 = [x0; Phi];                     %Initial conditions

% Generate initial trajectory 
tspan = 0:dt:T; 
[~, x1] = ode113(@(t,x)volterra_dynamics(t, x, param), tspan, x0, options);
s = x1;

% Differential correction scheme 
tol = 1e-5;         %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 10;       %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(2, maxIter);

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = (s(end,1:2)-s(1,1:2)).'; 
    
    %Evaluate the STM at the final state 
    STM = reshape(s(end,3:end), 2, 2); 
    
    %Compute the needed correction 
    dx0(:,iter) = (STM-eye(2))\e;
    
    %Generate new trajectory 
    x0(1:2) = x0(1:2)-dx0(1:2,iter);                                           %New initial conditions
    [~, s] = ode113(@(t,x)volterra_dynamics(t, x, param), tspan, x0, options); %Dynamics integration
    
    %Convergence 
    if (norm(e) < tol)
        GoOn = false;
    else
        iter = iter+1;
    end
end

% Results
if (~GoOn)
    disp('Case study 3: PASS.');
else
    disp('Case study 3: FAIL. Error convergence was not met.');
end
 
% Plotting
figure(3) 
hold on
plot(x1(:,1), x1(:,2), '-.r');
plot(s(:,1), s(:,2), 'b');
hold off 
xlabel('Position'); 
ylabel('Velocity');
title('Differential correction results in the phase space'); 
legend('Initial trayectory', 'Final trayectory');

%% Scenario 4: forcing von Karman street at desired frequency.
% Results: impossible solution. The frequency of the oscilation is
% unconstrained by the model. Only the amplitude of the oscillation could be
% corrected, but no initial conditions targetting is needed.

% Set up integration
options = odeset('RelTol', 2.25e-10, 'AbsTol', 1e-22);

% Define Strouhal number 
St = 1;
T = 1/St;

% Augmented initial conditions 
[dt, x0, param] = initial_conditions2(); 
Phi = eye(6);                               %Initial STM
Phi = reshape(Phi, [36 1]);
x0 = [x0; Phi];                             %Initial conditions

% Generate initial trajectory 
tspan = 0:dt:T; 
[~, x1] = ode113(@(t,x)karman_dynamics(t, x, param), tspan, x0, options);
s = x1;

% Differential correction scheme 
tol = 1e-5;         %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 1;        %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(6, maxIter);

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = (s(end,1:6)-s(1,1:6)).'; 
    
    %Evaluate the STM at the final state
    STM = reshape(s(end,7:end), 6, 6);      
    
    %Compute the needed correction 
    dx0(:,iter) = STM\e;
    
    %Generate new trajectory 
    x0(1:6) = x0(1:6)-dx0(:,iter);                                            %New initial conditions                                     
    [~, s] = ode113(@(t,x)karman_dynamics(t, x, param), tspan, x0, options);  %Dynamics integration
    
    %Convergence 
    if (norm(e) < tol)
        GoOn = false;
    else
        iter = iter+1;
    end
end

% Results
if (~GoOn)
    disp('Case study 4: PASS.');
else
    disp('Case study 4: FAIL. Error convergence was not met.');
end

%% Auxiliary functions
% Generate initial conditions for Duffing dynamics
function [dt, x, param, T] = initial_conditions()
    %Constants
    dt = 1e-3;  %Time step
    
    %Dynamical parameters 
    param(1) = 1;
    param(2) = 2; 
    param(3) = 0.02;
    param(4) = 1;
    param(5) = 0.5;
    
    %Select period randomly 
    T = 2*pi/param(5);
    
    %Select random initial conditions for both position and velocity 
    x = rand(2,1);
end

% Generate initial conditions for von Karman street ROM dynamics
function [dt, x, param] = initial_conditions2()
    %Constants
    dt = 1e-3;  %Time step
    
    %Dynamical parameters 
    param(1) = 1/10;
        
    %Select random initial conditions for both position and velocity 
    x = rand(6,1);
    x(3) = 0; 
    x(6) = 0;
end

% Generate initial conditions for the Lorenz dynamics
function [dt, x, param, T] = initial_conditions3()
    %Constants
    dt = 1e-3;  %Time step
    
    %Dynamical parameters 
    param(1) = 10;
    param(2) = 28;
    param(3) = 8/3;
    
    %Orbit period
    T = 6.02;
        
    %Select initial conditions 
    x = [0.695295; 1.27032; 14.132];
end

% Generate initial conditions for the Lokta-Volterra dynamics 
function [dt, x, param, T] = initial_conditions4()
    %Constants
    dt = 1e-3;  %Time step
    
    %Dynamical parameters 
    param = rand(4,1);
    
    %Orbit period
    T = 2*pi/(param(1)*param(3));
        
    %Select initial conditions 
    x = rand(2,1);
end

% Dynamics of the Duffing equation
function [ds] = duffing_dynamics(t, s, param)
    %Dynamical variables 
    x = s(1);                               %Position
    v = s(2);                               %Velocity
    Phi0 = reshape(s(3:end), [2 2]);        %State transition matrix
    
    %Jacobian of the system 
    J = [0 1; -param(1)-3*param(2)*x^2 -param(3)];
    
    %Variational equations
    Phi = J*Phi0;
    Phi = reshape(Phi, [4 1]);
    
    %Vector field 
    ds = [v; 
          -param(3)*v-param(1)*x-param(2)*x^3+param(4)*cos(param(5)*t); 
          Phi];
end

% Dynamics of the vortex shedding ROM
function [ds] = karman_dynamics(~, s, param)
    %Dynamical variables 
    u = s(1);                           %First principal mode
    v = s(2);                           %Second principal mode
    w = s(3);                           %Third principal mode    
    Phi = reshape(s(7:end), [6 6]);     %Initial STM
    mu = param(1);                      %Parameter of the system
    
    %Jacobian of the system
    H = [mu -1 -u; 1 mu -v; 2*u 2*v -1];
    J = [zeros(3,3) eye(3); H zeros(3,3)];
    
    %Variational equations 
    Phi = J*Phi; 
    Phi = reshape(Phi, [36 1]); 
    
    %Vector field 
    ds = [    s(4:6);
          mu*u-v-u*w; 
          u+mu*v-v*w; 
          u^2+v^2-w; 
                 Phi];    
end

% Dynamics of the Lorenz system
function [ds] = lorenz_dynamics(~, s, param)
    %Constants of the model 
    a = param(1);       %Prandtl number
    b = param(2);       %Rayleigh number
    c = param(3);
    
    %Variables of the problem 
    x = s(1); 
    y = s(2); 
    z = s(3);
    Phi = reshape(s(4:end), 3, 3);     %Initial STM 
    
    %Jacobian of the system
    J = [-a a 0; b -1 -x; y x -c];
    
    %Variational equations 
    Phi = J*Phi; 
    Phi = reshape(Phi, [9 1]); 
    
    %Vector field 
    ds = [a*(y-x); 
          x*(b-z)-y; 
          x*y-c*z;
               Phi]; 
end

% Dynamics of the Lotka-Volterra system 
function [ds] = volterra_dynamics(~, s, param)
    %Dynamical constants 
    n = 2;              %Phase space dimension
    alpha = param(1); 
    beta = param(2); 
    gamma = param(3); 
    delta = param(4); 
    
    %State variables
    x = s(1);       %Number of preys
    y = s(2);       %Number of predatos 
    
    %Variational equations
    Phi = reshape(s(n+1:end), n, n);        %Initial STM 
    J = [alpha -beta*x; gamma*y -delta];    %Jacobian of the system
    Phi = J*Phi;                            %Variational equations
    Phi = reshape(Phi, [n^2 1]);
    
    %Vector field
    ds = [x*(alpha-beta*y); 
         -y*(gamma-delta*x);
                       Phi];
    
end