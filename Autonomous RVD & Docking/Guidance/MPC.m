%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 6: MPC guidance-control law %% 
% This script provides an interface to test MPC control law rendezvous strategies for
% rendezvous missions.

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
%Phase space dimension 
n = 6; 

%Time span 
dt = 1e-3;                          %Time step
tf = 0.6;                           %Rendezvous time
tspan = 0:dt:tf;                    %Integration time span
tspann = 0:dt:2*pi;                 %Integration time span

%CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
L = libration_points(mu);           %System libration points
Lem = 384400e3;                     %Mean distance from the Earth to the Moon

%Differential corrector set up
nodes = 10;                         %Number of nodes for the multiple shooting corrector
maxIter = 20;                       %Maximum number of iterations
tol = 1e-10;                        %Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
%Halo characteristics 
Az = 200e6;                                                         %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                                     %Bifucartion tolerance on the stability index
num = 2;                                                            %Number of orbits to continuate
method = 'SPC';                                                     %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                        %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};                 %Object and characteristics to continuate
corrector = 'Plane Symmetric';                                      %Differential corrector method
direction = 1;                                                      %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                                 %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(2,:), maxIter, tol);

%% Modelling in the synodic frame %%
index = fix(tf/dt);                                         %Rendezvous point
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[~, S] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% Optimal control guidance scheme
%Set up of the optimization
method = 'Genetic algorithm';                 %Method to solve the problem
impulses = 5;                                 %Number of impulses
TOF = tspan(end);                             %Time of flight
cost_function = 'State';                   %Target a position rendezvous

%Thruster characteristics 
Tmin = -0.1;                                  %Minimum thrust capability (in velocity impulse)
Tmax = 0.1;                                   %Maximum thrust capability (in velocity impulse)

%Main computation 
[St, dV, state] = OPTI(mu, cost_function, Tmin, Tmax, TOF, s0, impulses, method);

dVl1(1:3,1) = sum(dV,2);                      %L1 norm of the impulses 
dVl2(1) = sum(sqrt(dot(dV,dV,2)));            %L2 norm of the impulses 

%Error in time 
e = zeros(1,size(St,1));                      %Preallocation of the error
for i = 1:size(St,1)
    e(i) = norm(St(i,7:12));
end
e(1) = norm(Sn(1,7:12));                      %Initial error before the burn

%Compute the error figures of merit 
ISE = trapz(tspan, e.^2);
IAE = trapz(tspan, abs(e));

%% Results %% 
disp('SIMULATION RESULTS: ')
if (state.State)
    disp('   Multi impulsive rendezvous was achieved');
    fprintf('   Delta V budget (L1 norm): %.4ei %.4ej %.4ek \n', dVl1(1,1), dVl1(2,1), dVl1(3,1));
    fprintf('   Delta V budget (L2 norm): %.4e \n', dVl2(:,1));
else
    disp('    Multi impulsive rendezvous was not achieved');
end

%Plot results 
figure(1) 
view(3) 
hold on
plot3(Sn(:,1), Sn(:,2), Sn(:,3)); 
plot3(S_rc(:,1), S_rc(:,2), S_rc(:,3)); 
hold off
legend('Target motion', 'Chaser motion'); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Reconstruction of the natural chaser motion');

%Plot relative phase trajectory
figure(2) 
view(3) 
plot3(St(:,7), St(:,8), St(:,9)); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative motion in the configuration space');

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan, St(:,7)); 
plot(tspan, St(:,8)); 
plot(tspan, St(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('x coordinate', 'y coordinate', 'z coordinate');
title('Relative position evolution');
subplot(1,2,2)
hold on
plot(tspan, St(:,10)); 
plot(tspan, St(:,11)); 
plot(tspan, St(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinate');
grid on;
legend('x velocity', 'y velocity', 'z velocity');
title('Relative velocity evolution');

%Configuration space error 
figure(4)
plot(tspan, log(e)); 
xlabel('Nondimensional epoch');
ylabel('Absolute error  (log)');
grid on;
title('Absolute error in the configuration space (L2 norm)');

%Rendezvous animation 
if (false)
    figure(5) 
    view(3) 
    grid on;
    hold on
    plot3(St(:,1), St(:,2), St(:,3), 'k-.'); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    title('Rendezvous simulation');
    for i = 1:size(St,1)
        T = scatter3(St(i,1), St(i,2), St(i,3), 30, 'b'); 
        V = scatter3(St(i,1)+St(i,7), St(i,2)+St(i,8), St(i,3)+St(i,9), 30, 'r');

        drawnow;
        delete(T); 
        delete(V);
    end
    hold off
end

%% Auxiliary functions
%Optimization function 
function [Sg, dV, state] = OPTI(mu, cost_function, Tmin, Tmax, TOF, s0, impulses, method)
    %Constants 
    m = 6;                                  %Phase space dimension
    
    %Sanity check on the dimension 
    if (size(s0,1) ~= 2*m)
        s0 = s0.';                          %Initial conditions
    end
    
    %Differential corrector setup
    tol = 1e-8;                             %Differential corrector setup 
    maxIter = 50;                           %Maximum number of iterations
    iter = 1;                               %Initial iteration
    GoOn = true;                            %Convergence boolean
    
    %Integration setup 
    RelTol = 2.25e-14; 
    AbsTol = 1e-22; 
    options = odeset('RelTol', RelTol, 'AbsTol', AbsTol);
    
    %Initial integration    
    Phi = eye(m);                           %Initial STM
    Phi = reshape(Phi, [m^2 1]);            %Reshape the initial STM
    s0 = [s0; Phi];                         %Complete initial conditions
    dt = 1e-3;                              %Time step
    tspan = 0:dt:TOF;                       %Integration time span
    
    [~, Sn] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), tspan, s0, options);
    St = Sn; 
    
    %Preallocation 
    dV = zeros(maxIter,3,length(tspan));    %Velocity impulses preallocation all along the TOF orbit arc
    
    %Main computation 
    while ((GoOn) && (iter < maxIter))
        %Compute the commands 
        [commands, time_indexes] = opt_core(cost_function, Tmin, Tmax, tspan(end), dt, St, impulses, method); 
        time_indexes = sort(fix(time_indexes/dt)+1);        %Sort the firings times 
        
        k = 1;                                              %Target initial maneuver
        for i = 1:length(tspan)
            if (k <= impulses)
                if (i == time_indexes(k))
                    dV(iter,:,i) = commands(:,k);           %Save the firing at each particular moment
                    k = k+1;                                %Target the new firing
                end
            end
        end
        
        %Recompute the trajectory 
        for i = 1:length(tspan)-1
            %Add the maneuver
            St(i,10:12) = St(i,10:12) + shiftdim(dV(iter,:,i)).';  
            
            %New integration
            [~, Saux] = ode113(@(t,s)nlr_model(mu, true, false, true, 'Encke', t, s), [0 dt], St(i,:), options);
            
            %Next initial conditions
            St(i+1,:) = Saux(end,:);
        end
        
        %Compute the rendezvous error
        switch (cost_function)
            case 'Position'
                e = St(end,7:9).';      
            case 'Velocity'
                e = St(end,10:12).'; 
            case 'State'
                e = St(end,7:12).'; 
            otherwise
                error('No valid cost function was chosen');
        end
        
        %Convergence analysis
        norm(e)
        if (norm(e) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end
    
    %Output 
    dV = shiftdim(sum(dV,1));       %Control law at each discrete point over time
    Sg = St;                        %Optimal control trajectory
    state.State = ~GoOn;            %Convergence state
    state.Error = norm(e);          %Final error 
    state.Iterations = iter;        %Final number of iterations
end

%Core nonlinear optimization function
function [commands, time_indexes] = opt_core(cost_function, Tmin, Tmax, TOF, dt, trajectory, impulses, method)    
    %Linear constraints 
    A = []; 
    b = []; 
    Aeq = []; 
    beq = [];
    
    %Upper and lower bounds
    lb = [zeros(1,impulses) Tmin*ones(1,3*impulses)];         %Lower bound
    ub = [TOF*ones(1,impulses) Tmax*ones(1,3*impulses)];      %Upper bound
    
    switch (method)
        case 'Genetic algorithm'
            %General set up
            dof = impulses*(3+1);   %Three-dimensional control vector for each instant
            PopSize = 100;         %Population size for each generation
            MaxGenerations = 10;   %Maximum number of generations for the evolutionary algorithm
            
            options = optimoptions(@ga,'PopulationSize', PopSize, 'MaxGenerations', MaxGenerations);
                            
            %Compute the commands
            solution = ga(@(u)costfunc(impulses,u), dof, A, b, Aeq, beq, lb, ub, ...
                          @(u)nonlcon(cost_function, impulses, dt, trajectory, u), options);
            
            time_indexes = solution(1,1:impulses);                            %Times at which the firings are perfomed
            commands = reshape(solution(1,impulses+1:end), 3, impulses);      %Control law
            
        case 'NPL'
            %Initial guess
            sol0 = [zeros(1,impulses) Tmax*ones(1,3*impulses)];        
            
            %Compute the commands
            solution = fmincon(@(u)costfunc(impulses,u), sol0, A, b, Aeq, beq, lb, ub, ...
                               @(u)nonlcon(cost_function, impulses, dt, trajectory, u));
            
            time_indexes = solution(1,1:impulses);                            %Times at which the firings are perfomed
            commands = reshape(solution(1,impulses+1:end), 3, impulses);      %Control law
            
        otherwise 
            error('No valid method was chosen');
    end
end

%Cost function 
function [cost] = costfunc(impulses,u)
    %Reshape the control law 
    u = reshape(u(1,impulses+1:end), [3 impulses]);
    
    %Cost function 
    cost = sum(sum(abs(u),1));
end

%Nonlinear constraints
function [c, ceq] = nonlcon(cost_function, impulses, dt, trajectory, x)
    %Constants 
    m = 6;                  %Phase space dimension
    tol = 1e-5;             %Rendezvous tolerance
    tol = tol*ones(m,1);    %State tolerance
    
    %Natural STM map
    Monodromy = reshape(trajectory(end,2*m+1:end), [m m]);          %Final STM
    error = trajectory(end,m+1:2*m).';                              %Final state, also the error to rendezvous
    
    %Index the impulsive times 
    times = fix(x(1,1:impulses)/dt)+1;
    times = sort(times);
    
    %Generate the inequality function
    control = zeros(m,1);                                           %Total control effort
    for i = 1:impulses
        index = times(i);                                           %Time index to perform the firings
        u = reshape(x(1,impulses+1+3*(i-1):impulses+3*i),3,1);      %Control law
        STM = reshape(trajectory(index,2*m+1:end), [m m]);          %STM from the initial time to the time ti
        STM = Monodromy*STM^(-1);                                   %Relative STM from time ti to TOF
        control = control + STM(:,4:6)*u;                           %Accumulated control effort                       
    end
    
    %Nonlinear constraints
    c = control+error-tol;                  %Target the complete rendezvous
    switch (cost_function)
        case 'Position'
            c = c(1:3);                     %Target a position rendezvous
        case 'Velocity'
            c = c(4:6);                     %Target a velocity rendezvous
        case 'State'

        otherwise 
            error('No valid cost function was chosen');
    end
    
    ceq = [];                               %Empty equality constraint
end

%Core linear optimization function
function [commands] = lopt_core(cost_function, Tmin, Tmax, trajectory) 
    %Constants 
    m = 6;                                                            %Phase space dimension
    dim = size(trajectory,1);                                         %Number of impulses to make 
    Monodromy = reshape(trajectory(end,2*m+1:end), [m m]);            %Final STM
    
    %Linear equality constraints 
    beq = trajectory(end,m+1:2*m).';                                  %Error to rendezvous
    
    switch (cost_function)
        case 'Position'
            Aeq = zeros(m/2, m*dim);                                  %Preallocate the equality linear constraint matrix
            for i = 1:dim
                STM = reshape(trajectory(i,2*m+1:end), [m m]);        %STM at each point
                STM = Monodromy*STM^(-1);                             %Relative STM
                Aeq(:,1+(m/2)*(i-1):(m/2)*i) = STM(1:3,4:6);          %Final equality linear constraint matrix   
            end
            beq = beq(1:3);                                           %Target position
            
        case 'Velocity'
            Aeq = zeros(m/2, m*dim);                                  %Preallocate the equality linear constraint matrix
            for i = 1:dim
                STM = reshape(trajectory(i,2*m+1:end), [m m]);        %STM at each point
                STM = Monodromy*STM^(-1);                             %Relative STM
                Aeq(:,1+(m/2)*(i-1):(m/2)*i) = STM(4:6,4:6);          %Final equality linear constraint matrix   
            end
            beq = beq(4:6);                                           %Target velocity
            
        case 'State'
            Aeq = zeros(m, m*dim);                                    %Preallocate the equality linear constraint matrix
            for i = 1:dim
                STM = reshape(trajectory(i,2*m+1:end), [m m]);        %STM at each point
                STM = Monodromy*STM^(-1);                             %Relative STM
                Aeq(:,1+(m/2)*(i-1):(m/2)*i) = STM(:,4:6);            %Final equality linear constraint matrix   
            end
            
        otherwise
            error('No valid cost function was chosen');
    end
        
    %Linear inequality constraints
    A = [eye(3*dim) zeros(3*dim);  zeros(3*dim), -eye(3*dim)];
    A(1:3*dim,3*dim+1:end) = -eye(3*dim);
    A(3*dim+1:end,1:3*dim) = -eye(3*dim);
    b = zeros(m*dim,1);
    
    %Upper and lower bounds
    lb = [Tmin*ones(3*dim,1); zeros(3*dim,1)];             %Lower bound
    ub = [Tmax*ones(3*dim,1); abs(Tmax)*ones(3*dim,1)];    %Upper bound
    
    %Cost function 
    f = [zeros(1,3*dim) ones(1,3*dim)];
    
    %Solve the problem 
    sol = linprog(f,A,b,Aeq,beq,lb,ub);
    
    %Output 
    commands = reshape(sol(1:3*dim,1), [3 dim]);
end