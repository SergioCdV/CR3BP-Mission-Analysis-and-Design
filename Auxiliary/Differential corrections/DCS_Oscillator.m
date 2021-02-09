%% Differential correction example: simple linear oscillator %% 
% Sergio Cuevas del Valle % 
% December 2020 % 

%% Description %%
% This file provides several examples of differential corrections schemes applied to a simple 
% linear oscillator. Initial conditions, chosen randomly, are modified to
% accomplish different final conditions in both amplitude and wave speed.
% The file iteratively solves several 2 boundary problems by means of the
% state transition matrix of the system, known analytically in this case.

%% Scenario 1: target final position and velocity varying the whole state vector. 
% Results: convergence.

% Initial conditions 
[dt, tf, omega, x0, xf, phase] = initial_conditions(); 

% Generate initial trajectory 
t = 0:dt:tf; 
x1 = [x0(1)*cos(omega*t+phase(1))+(x0(2)/omega)*sin(omega*t+phase(2));    %Initial position trajectory
     -x0(1)*omega*sin(omega*t+phase(1))+x0(2)*cos(omega*t+phase(2))];     %Initial velocity trajectory
x = x1;

% Differential correction scheme 
tol = 1e-5;         %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 1e3;      %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(2, maxIter);

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = x(:,end)-xf; 
    
    %Evaluate the STM at the final state 
    STM = [cos(omega*tf+phase(1)) (1/omega)*sin(omega*tf+phase(2)); 
          -omega*sin(omega*tf+phase(1)) cos(omega*tf+phase(2))]; 
    
    %Compute the needed correction 
    dx0(:,iter) = STM\e;
    
    %Generate new trajectory 
    x0 = x0-dx0(:,iter);                                                     %New initial conditions
    x = [x0(1)*cos(omega*t+phase(1))+(x0(2)/omega)*sin(omega*t+phase(2));    %Position trajectory
        -x0(1)*omega*sin(omega*t+phase(1))+x0(2)*cos(omega*t+phase(2))];     %Velocity trajectory
    
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
plot(t, x1(1,:), '-.r');
plot(t, x(1,:), 'b');
plot(tf, xf(1), 'o');
hold off 
xlabel('Time (s)'); 
ylabel('Amplitude (s)');
title('Differential correction results for an oscillator'); 
legend('Initial trayectory', 'Final trayectory', 'Desired final amplitude');

%% Scenario 2: target final position and velocity varying only the initial position vector.
% Results: no convergence. Not enough degrees of freedom.

% Initial conditions 
[dt, tf, omega, x0, xf, phase] = initial_conditions();

% Generate initial trajectory 
t = 0:dt:tf; 
x = [x0(1)*cos(omega*t+phase(1))+(x0(2)/omega)*sin(omega*t+phase(2));    %Initial position trajectory
    -x0(1)*omega*sin(omega*t+phase(1))+x0(2)*cos(omega*t+phase(2))];     %Initial velocity trajectory

% Differential correction scheme 
tol = 1e-5;         %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 1e3;      %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(1, maxIter);

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = x(:,end)-xf; 
    
    %Evaluate the STM at the final state 
    STM = [cos(omega*tf+phase(1)) (1/omega)*sin(omega*tf+phase(2)); 
          -omega*sin(omega*tf+phase(1)) cos(omega*tf+phase(2))]; 
    
    %Compute the needed correction 
    dx0(iter) = pinv(STM(:,1))*e;
    
    %Generate new trajectory 
    x0(1) = x0(1)-dx0(iter);                                                 %New initial conditions
    x = [x0(1)*cos(omega*t+phase(1))+(x0(2)/omega)*sin(omega*t+phase(2));    %Position trajectory
        -x0(1)*omega*sin(omega*t+phase(1))+x0(2)*cos(omega*t+phase(2))];     %Velocity trajectory
    
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

%% Scenario 3: target final position and velocity varying only the initial velocity vector.
% Results: no convergence. Not enough degrees of freedom.

% Initial conditions 
[dt, tf, omega, x0, xf, phase] = initial_conditions();

% Generate initial trajectory 
t = 0:dt:tf; 
x = [x0(1)*cos(omega*t+phase(1))+(x0(2)/omega)*sin(omega*t+phase(2));    %Initial position trajectory
    -x0(1)*omega*sin(omega*t+phase(1))+x0(2)*cos(omega*t+phase(2))];     %Initial velocity trajectory

% Differential correction scheme 
tol = 1e-5;         %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 1e3;      %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(1, maxIter);

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = x(:,end)-xf; 
    
    %Evaluate the STM at the final state 
    STM = [cos(omega*tf+phase(1)) (1/omega)*sin(omega*tf+phase(2)); 
          -omega*sin(omega*tf+phase(1)) cos(omega*tf+phase(2))]; 
    
    %Compute the needed correction 
    dx0(iter) = pinv(STM(:,2))*e;
    
    %Generate new trajectory 
    x0(2) = x0(2)-dx0(iter);                                                 %New initial conditions
    x = [x0(1)*cos(omega*t+phase(1))+(x0(2)/omega)*sin(omega*t+phase(2));    %Position trajectory
        -x0(1)*omega*sin(omega*t+phase(1))+x0(2)*cos(omega*t+phase(2))];     %Velocity trajectory
    
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

%% Scenario 4: varying the integration time and position to obtain a peak
% Results: general convergence.

% Initial conditions 
[dt, ~, omega, x0, ~, phase] = initial_conditions();

% Target final conditions
xf = [sqrt(x0(1)^2+(x0(2)/omega)^2+2*cos(phase(2))*(x0(2)/omega)*x0(1)); 0];

% Generate initial trajectory for plotting 
tf = (1/5)*(2*pi/omega);
t = 0:dt:tf; 
x = [x0(1)*cos(omega*t+phase(1))+(x0(2)/omega)*sin(omega*t+phase(2));    %Initial position trajectory
    -x0(1)*omega*sin(omega*t+phase(1))+x0(2)*cos(omega*t+phase(2))];     %Initial velocity trajectory

% Differential correction scheme 
tol = 1e-4;         %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 1e3;      %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(2, maxIter);

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = x(:,end)-xf; 
    
    %Evaluate the vector field at tf
    F = [x(2,end); -omega^2*x(1,end)]; 
    
    %Compute the augmented STM 
    STM = [ cos(omega*tf+phase(1))      F(1); 
          -omega*sin(omega*tf+phase(1)) F(2)];
          
    %Compute the needed correction 
    dx0(:,iter) = STM\e;    
    
    %Generate new trajectory 
    x0(1) = x0(1)-dx0(1,iter);                                               %New initial position
    tf = tf-dx0(2,iter);                                                     %New integration time
    t = 0:dt:tf;                                                             %New time span
    x = [x0(1)*cos(omega*t+phase(1))+(x0(2)/omega)*sin(omega*t+phase(2));    %Position trajectory
        -x0(1)*omega*sin(omega*t+phase(1))+x0(2)*cos(omega*t+phase(2))];     %Velocity trajectory
    
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

%% Scenario 5: varying the integration time and the velocity vector to obtain a peak
% Results: general convergence.

% Initial conditions 
[dt, ~, omega, x0, ~, phase] = initial_conditions();

%Target final conditions
xf = [sqrt(x0(1)^2+(x0(2)/omega)^2 +2*cos(phase(2))*(x0(2)/omega)*x0(1)); 0];

% Generate initial trajectory
tf = (1/5)*(2*pi/omega);
t = 0:dt:tf; 
x = [x0(1)*cos(omega*t+phase(1))+(x0(2)/omega)*sin(omega*t+phase(2));    %Initial position trajectory
    -x0(1)*omega*sin(omega*t+phase(1))+x0(2)*cos(omega*t+phase(2))];     %Initial velocity trajectory

% Differential correction scheme 
tol = 1e-4;         %Tolerance to stop the process
iter = 1;           %Initial iteration
maxIter = 1e3;      %Maximum number of iterations
GoOn = true;        %Convergence flag 

% Preallocation
dx0 = zeros(2, maxIter);

% Main computation
while (GoOn) && (iter < maxIter)
    %Compute the error vector 
    e = x(:,end)-xf; 
    
    %Evaluate the vector field at tf
    F = [x(2,end); -omega^2*x(1,end)]; 
    
    %Evaluate the STM at the final state 
    STM = [(1/omega)*sin(omega*tf+phase(2)) F(1); 
           cos(omega*tf+phase(2))           F(2)];  
      
    %Compute the needed correction 
    dx0(:,iter) = STM\e;
    
    %Generate new trajectory 
    x0(2) = x0(2)-dx0(1,iter);                                               %New initial conditions
    tf = tf-dx0(2,iter);                                                     %New integration time
    t = 0:dt:tf;                                                             %New time span
    x = [x0(1)*cos(omega*t+phase(1))+(x0(2)/omega)*sin(omega*t+phase(2));    %Position trajectory
        -x0(1)*omega*sin(omega*t+phase(1))+x0(2)*cos(omega*t+phase(2))];     %Velocity trajectory
    
    %Convergence 
    if (norm(e) < tol)
        GoOn = false;
    else
        iter = iter+1;
    end
end

% Results
if (~GoOn)
    disp('Case study 5: PASS.');
else
    disp('Case study 5: FAIL. Error convergence was not met.');
end

%% Auxiliary functions 
function [dt, tf, omega, x0, xf, phase] = initial_conditions()
    dt = 1e-3;              %Time step
    tf = 10;                %Integration time in normalized coordinates 
    omega = rand;           %Random angular velocity
    x0 = rand(2,1);         %Initial conditions randomly chosen in normalized coordinates 
    xf = rand(2,1);         %Final conditions randomly chosen in normalized coordinates
    phase = [0; 0];         %Initial phase randomly chosen
end