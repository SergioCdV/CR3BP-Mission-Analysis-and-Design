%% Differential correction example: lunar mission return trayectory %% 
% Sergio Cuevas del Valle % 
% January 2021 % 

%% Description %%
% This file provides a differential correction scheme to generate 
% a free return trajectory for a Earth-Moon mission. Motion is confined 
% to the synodic plane, an CR3BP dynamics are used. Multiple 
% shooting differential corrections are used.

% Results: not working properly. Initial seed is correct but the differential correction schemes 
% fails to achieve the requirements of both Earth and Moon orbits.
% Reference to Pavlak Master Thesis, 2010

%% Constants and set up
% Initial conditions (to be input by the user) 
he = 200e3;                 %Parking orbit altitude at perigee
hl = 169e3;                 %Periselenum altitude at periselenum
alpha = 0;                  %Constraint on the initial fligth path angle

% Constants (normalizing magnitudes)
mu = 0.0121505;             %Reduced gravitational parameter of the system
L = 384403e3;               %Earth-Moon barycentric distance 
muE = 3.986e14;             %Gravitational parameter of the Earth
muM = 4.9048695e12;         %Moon gravitational parameter
m = 4;                      %Phase space dimension
n = 100;                    %Maximum number of iterations
Rsoi = L*(muM/muE)^(2/5);   %Lunar SOI radius  
T = 1/sqrt((muE+muM)/L^3);  %Earth-Moon synodic period
V = L/T;                    %Normalizing velocity of the Earth-Moon system

% Parking orbit elements
rpe = 6378e3+he;            %Initial parking orbit radius 

% Circumlunar orbit elements
rpl = 1738e3+hl;            %Periselenum altitude 

% Initial anomaly on the parking orbit at the departure. 
% Constrained dynamically by the synodic frame rotation at SOI crossing epoch
gamma = deg2rad(270)*ones(1,3);

%% Integration setup 
% Time span
dt = 1e-3;                  %Time step
tspan = 0:dt:100;           %Initial time span 

% Set up integration
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-14, 'Events', @(t,x)x_crossing(t,x));

%% Bisection method to find the SOI entry angle to produce a free return trajectory
% Set up method 
tol = 1e-10; 

% Initial guess and Earth ellipse transfer orbit elements
beta(1) = deg2rad(180);                           %SOI entry angle maximum value 
beta(2) = deg2rad(0);                             %SOI entry angle minimum value               
r1(1) = sqrt(Rsoi^2+L^2-2*Rsoi*L*cos(beta(1)));   %Earth ellipse radius at the SOI crossing 
r1(2) = sqrt(Rsoi^2+L^2-2*Rsoi*L*cos(beta(2)));   %Earth ellipse radius at the SOI crossing 
v0(1) = sqrt(2*muE*(1/rpe-1/(rpe+r1(1))));        %Initial orbital velocity maximum value
v0(2) = sqrt(2*muE*(1/rpe-1/(rpe+r1(2))));        %Initial orbital velocity minimum value

theta(1) = gamma(1)+asin(sin(beta(1))*(Rsoi/r1(1)));
theta(2) = gamma(2)+asin(sin(beta(2))*(Rsoi/r1(2)));

% Transfer orbit elements in the synodic frame
R0 = [-mu+rpe/L*cos(theta(1)) -mu+rpe/L*cos(theta(2)); 
          rpe/L*sin(theta(1))     rpe/L*sin(theta(2))];    %Initial position in the CR3BP
V0 = [-v0(1)/V*sin(theta(1)) -v0(2)/V*sin(theta(2)); 
       v0(1)/V*cos(theta(1))  v0(2)/V*cos(theta(2))];      %Initial velocity in the CR3BP
Phi = [eye(m) eye(m)];                                     %Initial STM 
Phi = reshape(Phi, [m^2 2]);                               %Initial STM
S0 = [R0; V0; Phi];                                        %Initial state in the CR3BP 

% Initial integration 
[~, Smax] = ode113(@(t,x)dynamics(mu, t, x), tspan, S0(:,1), options);
[~, Smin] = ode113(@(t,x)dynamics(mu, t, x), tspan, S0(:,2), options);

% Bisection method 
GoOn = true;        %Convergence flag
iter = 1;           %Initial iteration
iterMax = 100;      %Maximum number of iterations

while (GoOn) && (iter < iterMax) && (abs(beta(1)-beta(2))/2 > tol)
    % New interval
    beta(3) = sum(beta(1:2))/2; 
    
    % New transfer orbit elements
    r1(1) = sqrt(Rsoi^2+L^2-2*Rsoi*L*cos(beta(1)));   %Earth ellipse radius at the SOI crossing 
    r1(2) = sqrt(Rsoi^2+L^2-2*Rsoi*L*cos(beta(2)));   %Earth ellipse radius at the SOI crossing 
    r1(3) = sqrt(Rsoi^2+L^2-2*Rsoi*L*cos(beta(3)));   %Earth ellipse radius at the SOI crossing
    v0(1) = sqrt(2*muE*(1/rpe-1/(rpe+r1(1))));        %Initial orbital velocity maximum value
    v0(2) = sqrt(2*muE*(1/rpe-1/(rpe+r1(2))));        %Initial orbital velocity minimum value
    v0(3) = sqrt(2*muE*(1/rpe-1/(rpe+r1(3))));        %Initial orbital velocity minimum value
    
    theta(1) = gamma(1)+asin(sin(beta(1))*(Rsoi/r1(1)));
    theta(2) = gamma(2)+asin(sin(beta(2))*(Rsoi/r1(2)));
    theta(3) = gamma(3)+asin(sin(beta(3))*(Rsoi/r1(3)));

    % New transfer orbit elements in the synodic frame 
    R0 = [-mu+rpe/L*cos(theta(1)) -mu+rpe/L*cos(theta(2)) -mu+rpe/L*cos(theta(3)); 
              rpe/L*sin(theta(1))     rpe/L*sin(theta(2))     rpe/L*sin(theta(3))];  %Initial position in the CR3BP
    V0 = [-v0(1)/V*sin(theta(1)) -v0(2)/V*sin(theta(2)) -v0(3)/V*sin(theta(3)); 
           v0(1)/V*cos(theta(1))  v0(2)/V*cos(theta(2))  v0(3)/V*cos(theta(3))];     %Initial velocity in the CR3BP
    Phi = [eye(m) eye(m) eye(m)];                                                    %Initial STM 
    Phi = reshape(Phi, [m^2 3]);                                                     %Initial STM
    S0 = [R0; V0; Phi];                                                              %Initial state in the CR3BP

    % Trajectory integration 
    [~, Smax] = ode113(@(t,x)dynamics(mu, t, x), tspan, S0(:,1), options);
    [~, Smin] = ode113(@(t,x)dynamics(mu, t, x), tspan, S0(:,2), options);
    [t, Snew] = ode113(@(t,x)dynamics(mu, t, x), tspan, S0(:,3), options);
        
    % Convergence analysis 
    if (abs(Snew(end,3)) == 0)
        GoOn = false;
    else
        % New root interval
        if (sign(Smin(end,3)) == sign(Snew(end,3)))
            beta(2) = beta(3);
        else
            beta(1) = beta(3);
        end
        
        %Update iterations 
        iter = iter+1;
    end
end

% Display results 
if (GoOn)
    disp('Bisection algorithm convergence was not achieved.');
else
    disp('Bisection algorithm convergence was achieved.');
end

%% Differential correction scheme
TF = t(end);        %Semiperiod of the orbit
seed = Snew.';      %Initial conditions for the differential corrector
[S, state] = differential_corrector(mu, seed, 50, tol, TF, rpe/L, rpl/L, alpha); 

%% Results 
% Plotting
figure(1)
hold on
plot(Snew(:,1), Snew(:,2), '.-r');
plot(S.Trajectory(:,1), S.Trajectory(:,2), 'b');
hold off
title('Apollo program type mission trajectory'); 
xlabel('Normalized x coordinate');
ylabel('Normalized y coordinate');
legend('Initial converged guess', 'Final trajectory');
title('Mission trajectory');

%% Auxiliary function 
% Orbit event: x-axis crossing
function [Pos, isterminal, dir] = x_crossing(~, x)
    %Event defining
    Pos = x(2);
    isterminal = 1; 
    dir = -1;           %Direction of the crossing
end

% C3RBP dynamics 
function [dr] = dynamics(mu, ~, s)
    %Constants 
    n = 4;                  %Phase space dimension 
    
    %Define the initial phase space vector
    x = s(1);               %Synodyc x coordinate
    y = s(2);               %Synodyc y coordinate  
    V = s(3:4);             %Synodic velocity vector
    
    %Relevant system parameters
    mu1 = 1-mu;             %First primary normalized position
    mu2 = mu;               %Second primary normalized position
    r1 = [(x+mu2); y];      %Relative position vector to the first primary
    r2 = [(x-mu1); y];      %Relative position vector to the secondary primary
    R1 = norm(r1);          %Distance to the first primary
    R2 = norm(r2);          %Distance to the secondary primary
    
    %Compute the time flow of the system 
    inAcc = [x+2*V(2); y-2*V(1)];                    %Inertial acceleration
    F = [V; inAcc-mu1/R1^3*r1-mu2/R2^3*r2];          %Time flow of the system
    
    %Compute the initial STM
    phi = reshape(s(n+1:end), [n n]);

    %First variations of the augmented potential function (Hessian of the potential)
    Ux = 1-(mu1/R1^3)*(1-3*((x+mu2)/R1)^2)-(mu2/R2^3)*(1-3*((x-mu1)/R2)^2); 
    Uy = 1-(mu1/R1^3)*(1-3*(y/R1)^2)-(mu2/R2^3)*(1-3*(y/R2)^2);
    Uxy = 3*y*((mu1/R1^5)*(x+mu2)+(mu2/R2^5)*(x-mu1));
    Uyx = Uxy;

    %Compute the first variational equations evaluated at the reference
    O = zeros(2,2); 
    I = eye(2);
    G = [Ux Uxy; Uyx Uy];
    K = [0 2; -2 0];
    Jacob = [O I; G K];                             %Jacobian of the system   
    dphi = Jacob*phi;                               %Variational equations
    dphi = reshape(dphi, [n^2 1]); 

    %Update the differential configuration space vector
    dr = [F; dphi];
end

% Differential corrector to generate the free-return trajectory
function [xf, state] = differential_corrector(mu, seed, n, tol, T, rpe, rpl, alpha)    
    %Constants 
    m = 4;                          %Phase space dimension 
    Phi = eye(m);                   %Initial STM  
    Phi = reshape(Phi, [m^2 1]);    %Initial STM 
    dt = 1e-5;                      %Integration time step
    nodes = 5;                     %Number of relevant nodes
    Dt = T/nodes;                   %Time step
    constraints = 4;                %Additional constraints to continuity of the trajectory
    
    %Prepare initial conditions
    internalSeed = zeros((m+1)*nodes-1,1);        %Preallocate internal patch points seeds 
    
    %Divide the orbit into the internal nodes
    h = fix(size(seed,2)/nodes)-1;
    for i = 1:nodes
        internalSeed(m*(i-1)+1:m*i) = seed(1:m,(i-1)*h+1);
        if (i ~= nodes)
            internalSeed(end-(nodes-1)+i) = Dt;
        end
    end    
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', @(t,x)x_crossing(t,x));            

    %Set up differential correction scheme
    GoOn = true;                                          %Convergence flag
    maxIter = n;                                          %Maximum number of iterations   
    iter = 1;                                             %Initial iteration
    
    %Preallocation 
    ds0 = zeros(size(internalSeed,1),maxIter);            %Vector containing the initial conditions correction
    e = zeros(m*(nodes-1)+constraints,1);                 %Continuity error vector  
    A = zeros(m*(nodes-1)+constraints, m*nodes);          %STM matrix
    B = zeros(m*(nodes-1)+constraints, nodes-1);          %Dynamics matrix
        
    %Main computation 
    while (GoOn) && (iter < maxIter)        
        for i = 1:nodes
            %Proceed with the integration
            if (i ~= nodes)
                tspan = 0:dt:internalSeed(end-(nodes-1)+i);  
            else
                tspan = 0:dt:Dt;
            end
            S0 = [shiftdim(internalSeed(m*(i-1)+1:m*i)); Phi];                   %Subarc initial conditions
            [~, S] = ode113(@(t,s)dynamics(mu, t, s), tspan, S0, options);       %Integration
            F = dynamics(mu, 0, S(end,:).');                                     %Vector field 
            
            %Build the covariance matrix                                         %Vector field matrix
            if (i == nodes)
                %Constraints on the periselenum altitude and transversal crossing conditions
                STM = reshape(S(end,m+1:end),[m m]);                             %Subarc STM
                A(end-2,end-m+1:end) = 2*[S(end,1)+mu-1 S(end,2) 0 0];           %Periselenum altitude
                A(end-1:end,end-m+1:end) = STM(2:3,:);                           %Transversal x crossing constraint matrix
            else
                %Continuity constraint
                A(m*(i-1)+1:m*i,m*(i-1)+1:m*i) = reshape(S(end,m+1:end),[m m]);  %Subarc STM
                A(m*(i-1)+1:m*i,m*i+1:m*(i+1)) = -eye(m);                        %Propagation STM between arcs 
                B(m*(i-1)+1:m*i,i) = F(1:m);                                     %Dynamics matrix
                if (i == 1)
                    %Constraints on perigee altitude and fligth path angle
                    A(end-4,1:m) = 2*[S(1,1)+mu S(1,2) 0 0];                     %Perigee altitude constraint
                    A(end-3,1:m) = [S(1,3:4)+mu S(1,1)+mu S(1,2)];               %Perigee altitude constraint
                end
            end     
            
            %Compute the error 
            if (i == nodes)
                e(end-2) = norm([S(end,1)+mu-1 S(end,2)])^2-rpl^2;               %Periselenum constraint
                e(end-1:end) = S(end,2:3);                                       %Transversal crossing constraint
            else
                e(m*(i-1)+1:m*i) = shiftdim(S(end,1:m).'-internalSeed(m*i+1:m*(i+1)));
                if (i == 1)
                    e(end-4) = norm([S(1,1)+mu S(1,2)])^2-rpe^2;                 %Perigee altitude constraint
                    e(end-3) = dot([S(1,1)+mu S(1,2)], S(1,3:4))-sin(alpha);     %Flight path angle
                end
            end
        end
        
        %Full covariance matrix 
        C = [A B];
                
        %Compute the correction in the under-determined case
        ds0(:,iter) = C.'*(C*C.')^(-1)*e;     
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            internalSeed = internalSeed-ds0(:,iter);    %Update initial conditions
            iter = iter+1;                              %Update iteration
        end       
    end
    
    %Integrate the whole trayectory
    tspan = 0:dt:2*pi;    
    seed = [shiftdim(internalSeed(1:m)); Phi];                  
    [t, S] = ode113(@(t,s)dynamics(mu, t, s), tspan, seed, options);
    backS = S(end:-1:1,:); 
    backS(:,2) = -backS(:,2);
    backS(:,4) = -backS(:,4);
    
    %Ouput corrected trajectory 
    xf.Trajectory = [S; backS];     %Trayectory
    xf.Period = t(end);             %Orbit period
    
    %Ouput differential correction scheme convergence results
    state = ~GoOn;
end