%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 2: Two-impulsive rendezvous via STM %% 
% This script provides an interface to test the two-impusilve rendezvous strategy using the STM of 
% the dynamical model. 

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% The first impulse is known as targetting, aimed to nullify the relative
% position between chaser and target after some flight time tf. The second
% impulse nullifies the relative velocity between the two to complete the
% rendezvous.

% In the relative phase space, the relative particle is driven to the
% origin of the synodic frame.

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Set up %%
%Set up graphics 
set_graphics();

%Integration tolerances (ode113)
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);  

%% Contants and initial data %% 
%Time span 
dt = 1e-3;                          %Time step
tf = 0.6;                           %Rendezvous time
tspann = 0:dt:2*pi;                 %Integration time span

%CR3BP constants 
mu = 0.0121505;                     %Earth-Moon reduced gravitational parameter
L = libration_points(mu);           %System libration points
Lem = 384400e3;                     %Mean distance from the Earth to the Moon

%Differential corrector set up
maxIter = 50;                       %Maximum number of iterations
tol = 1e-10;                        %Differential corrector tolerance

%% Initial conditions and halo orbit computation %%
%Halo characteristics 
Az = 200e6;                                                 %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 1;                                                     %Orbits around L1
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                             %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');  %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%% Modelling in the synodic frame %% 
index = fix(tf/dt);                                         %Rendezvous point
r_t0 = target_orbit.Trajectory(100,1:6);                    %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
s0 = [r_t0 rho0];                                           %Initial conditions of the target and the relative state

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% GNC: two impulsive rendezvous, multiple shooting scheme %%
%Differential corrector set up
S = S(1:index,:);                           %Restrict the time integration span
T = index*dt;                               %Flight time along the arc
nodes = 3;                                  %Number of nodes to compute
GoOn = true;                                %Convergence boolean 
iter = 1;                                   %Initial iteration 

cost = 'Position';                          %Targeting rendezvous

%Preallocation 
dV = zeros(3,maxIter);                      %Targeting impulse

%Implementation 
[S, state] = MS_rendezvous(mu, S, T, nodes, maxIter, tol, cost);      %Trajectory optimization
St = S.Trajectory;                                                    %Final computed trajectory

%% Results %% 
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
plot3(St(:,43), St(:,44), St(:,45)); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative motion in the configuration space');

%Rendezvous animation 
if (false)
    figure(3) 
    view(3) 
    grid on;
    hold on
    plot3(Sn(1:index,1), Sn(1:index,2), Sn(1:index,3), 'k-.'); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    title('Rendezvous simulation');
    for i = 1:size(St,1)
        T = scatter3(Sn(i,1), Sn(i,2), Sn(i,3), 30, 'b'); 
        V = scatter3(St(i,1)+St(i,7), St(i,2)+St(i,8), St(i,3)+St(i,9), 30, 'r');

        drawnow;
        delete(T); 
        delete(V);
    end
    hold off
end

%% Auxiliary functions %% 
function [xf, state] = MS_rendezvous(mu, seed, T, nodes, maxIter, tol, cost)
    %Constants 
    m = 12;                         %Complete phase space dimension 
    n = 6;                          %Individual phase space dimension
    
    %Sanity check on the initial seed dimensions
    if (size(seed,2) == m) || (size(seed,1) == m)
        if (size(seed,2) == m)
            seed = seed.';          %Accomodate new format
        end
    else
        disp('No valid initial conditions');
        xf = []; 
        state = false; 
        return;
    end
      
    %Constants 
    Phi = eye(n);                           %Initial STM  
    Phi = reshape(Phi, [n^2 1]);            %Initial STM 
    dt = 1e-4;                              %Integration time step
    h = fix(size(seed,2)/nodes)-1;          %Temporal index step
    Dt = T/nodes;                           %Time step between arcs
    
    switch (cost)
        case 'Position'
            constraints = 3;                %Additional constraints to continuity (targeting rendezvous)
        case 'Velocity'
            constraints = 3;                %Additional constraints to continuity (relative position rendezvous)
        otherwise
            constraints = 6;                %Additional constraints to continuity (full rendezvous)
    end
        
    %Preallocate internal patch points seeds 
    internalSeed = zeros((m+1)*nodes-1,1);        
    
    %Divide the orbit into the internal nodes
    for i = 1:nodes
        internalSeed(m*(i-1)+1:m*i) = seed(1:m,(i-1)*h+1);
        if (i ~= nodes)
            internalSeed(end-(nodes-1)+i) = Dt;
        end
    end    
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Integration conditions and tolerances                     
    direction = 1;                                              %Forward integration
    
    %Set up differential correction scheme
    GoOn = true;                                                %Convergence flag
    iter = 1;                                                   %Initial iteration
    
    %Preallocation 
    ds0 = zeros(size(internalSeed,1),maxIter);                  %Vector containing the initial conditions correction
    e = zeros(m*(nodes-1)+constraints,1);                       %Error vector  
    A = zeros(m*(nodes-1)+constraints, m*nodes);                %STM matrix
    B = zeros(m*(nodes-1)+constraints, nodes-1);                %Dynamics matrix
        
    %Main computation 
    while (GoOn) && (iter < maxIter)        
        for i = 1:nodes
            %Proceed with the integration
            if (i ~= nodes)
                tspan = 0:dt:internalSeed(end-(nodes-1)+i);  
            else
                tspan = 0:dt:Dt;
            end          
            S0 = shiftdim(internalSeed(m*(i-1)+1:m*i));
            S0 = [S0(1:n); Phi; S0(n+1:end); Phi];
            [~, S] = ode113(@(t,s)nlr_model(mu, direction, true, 'Encke V', t, s), tspan, S0, options);   %New trajectory
            F = nlr_model(mu, direction, false, 'Encke', 0, [S(end,1:n) S(end,n+n^2+1:2*n+n^2)].');       %Vector field
            
            %Build the covariance matrix                                       
            if (i ~= nodes)
                %Continuity constraint
                STM = [reshape(S(end,n+1:n+n^2),[n n]) zeros(n,n); ...
                       zeros(n,n) reshape(S(end,2*n+n^2+1:end),[n n])];                 %Subarc STM
                A(m*(i-1)+1:m*i,m*(i-1)+1:m*i) = STM;                                   %Subarc STM
                A(m*(i-1)+1:m*i,m*i+1:m*(i+1)) = -eye(m);                               %Continuity constraint matrix
                B(m*(i-1)+1:m*i,i) = F(1:end);                                          %Dynamics matrix
                
                %Compute the continuity error
                e(m*(i-1)+1:m*i) = shiftdim([S(end,1:n) S(end,n+n^2+1:2*n+n^2)].'-internalSeed(m*i+1:m*(i+1)));  
            else
                %Rendezvous error
                STM = reshape(S(end,2*n+n^2+1:end),[n n]);                               %Corresponding STM
                switch (cost)
                    case 'Position'
                        e(end-constraints+1:end) = shiftdim(S(end,n+n^2+1:n+n^2+3)).';   %Relative phase space vector to the origin
                        A(end-constraints+1:end,end-constraints+1:end) = STM(1:3,4:6);   %Constraint matrix
                    case 'Velocity'
                        e(end-constraints+1:end) = shiftdim(S(end,n+n^2+4:2*n+n^2)).';   %Relative phase space vector to the origin
                        A(end-constraints+1:end,end-constraints+1:end) = STM(4:6,4:6);   %Constraint matrix
                    otherwise
                        e(end-n+1:end) = shiftdim(S(end,n+n^2+1:2*n+n^2)).';             %Relative phase space vector to the origin
                        A(end-constraints+1:end,end-2:end) = STM(:,4:6);                 %Constraint matrix
                end
            end     
        end
        
        %Full covariance matrix 
        C = [A B];
                
        %Compute the correction 
        ds0(:,iter) = C.'*(C*C.')^(-1)*e;               %Compute the variation (under-determined case)
        
        %Convergence analysis 
        norm(e)
        if (norm(e) <= tol)
            GoOn = false;
        else
            internalSeed = internalSeed-ds0(:,iter);    %Update initial conditions
            iter = iter+1;                              %Update iteration
        end       
    end
    
    %Integrate the whole trayectory
    tspan = 0:dt:sum(internalSeed(end-nodes+1:end))+Dt; 
    seed = [shiftdim(internalSeed(1:n)); Phi; shiftdim(internalSeed(n+1:2*n)); Phi];
    [t, S] = ode113(@(t,s)nlr_model(mu, direction, true, 'Encke V', t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;                           %Trajectory
    xf.Period = t(end);                          %Orbit period
    
    if (GoOn)
        xf.Impulses = ds0(:,iter);               %Needed impulses to correct the trajectory
    else
        xf.Impulses = ds0(:,iter-1);             %Needed impulses to correct the trajectory
    end
        
    %Ouput differential correction scheme convergence results
    state = ~GoOn;
end

