%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 3: LQR/SDRE control law %% 
% This script provides an interface to test the LQR rendezvous strategy for
% rendezvous missions

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% The classical LQR for a time-varying model is solved to drive the relative phase space vector to the origin 
% (rendezvous condition).

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

%Spacecraft mass 
mass = 1e-10;

%Time span 
dt = 1e-3;                          %Time step
tf = 0.4;                           %Rendezvous time
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

%Continuate the first halo orbit to locate the chaser spacecraft
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 2;                                                    %Number of orbits to continuate
method = 'SPC';                                             %Type of continuation method (Single-Parameter Continuation)
algorithm = {'Energy', NaN};                                %Type of SPC algorithm (on period or on energy)
object = {'Orbit', halo_seed, target_orbit.Period};         %Object and characteristics to continuate
corrector = 'Plane Symmetric';                              %Differential corrector method
direction = 1;                                              %Direction to continuate (to the Earth)
setup = [mu maxIter tol direction];                         %General setup

[chaser_seed, state_PA] = continuation(num, method, algorithm, object, corrector, setup);
[chaser_orbit, ~] = differential_correction('Plane Symmetric', mu, chaser_seed.Seeds(2,:), maxIter, tol);

%% Modelling in the synodic frame %%
index = fix(tf/dt);                                         %Rendezvous point
if (index > size(target_orbit.Trajectory,1))
    index = mod(index, size(target_orbit.Trajectory,1));    %Rendezvous point
end
r_t0 = target_orbit.Trajectory(1,1:6);                      %Initial target conditions
r_c0 = target_orbit.Trajectory(1,1:6);                      %Initial chaser conditions 
rho0 = r_c0-r_t0;                                           %Initial relative conditions
rho0(1) = rho0(1)+100e3/Lem;
s0 = [r_t0 rho0].';                                         %Initial conditions of the target and the relative state

%Integration of the model
[t, S] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke', t, s), tspann, s0, options);
Sn = S;                

%Reconstructed chaser motion 
S_rc = S(:,1:6)+S(:,7:12);                                  %Reconstructed chaser motion via Encke method

%% Controlability analysis
%Approximation 
order = 2;                                   %Order of the approximation 

%Preallocation 
controlable = zeros(size(S(1:index),1),1);   %Controllability boolean

%Model coefficients 
mup(1) = 1-mu;                      %Reduced gravitational parameter of the first primary 
mup(2) = mu;                        %Reduced gravitational parameter of the second primary 
R(:,1) = [-mu; 0; 0];               %Synodic position of the first primary
R(:,2) = [1-mu; 0; 0];              %Synodic position of the second primary

%Linear model matrices
B = [zeros(n/2); zeros(n/2); eye(n/2)];         %Linear model input matrix 
Omega = [0 2 0; -2 0 0; 0 0 0];                 %Coriolis dyadic

for i = 1:size(S(1:index,:),1)
    %State coefficients 
    r_t = S(i,1:3).';               %Synodic position of the target
        
    %Relative position between the primaries and the target 
    Ur1 = r_t-R(:,1);               %Position of the target with respect to the first primary
    ur1 = Ur1/norm(Ur1);            %Unit vector of the relative position of the target with respect to the primary
    Ur2 = r_t-R(:,2);               %Position of the target with respect to the first primary
    ur2 = Ur2/norm(Ur2);            %Unit vector of the relative position of the target with respect to the primary
    
    %Compute the relative Legendre coefficient c2 
    rc = relegendre_coefficients(mu, S(i,1:3).', order); 
    cn = legendre_coefficients(mu, Ln, gamma, order);
    c2 = cn(2); 
    
    %Evaluate the linear model 
    %Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];
    Sigma = -((mup(1)/norm(Ur1)^3)+(mup(2)/norm(Ur2))^3)*eye(3)+3*((mup(1)/norm(Ur1)^3)*(ur1*ur1.')+(mup(2)/norm(Ur2)^3)*(ur2*ur2.'));
    A = [zeros(3) eye(3) zeros(3); zeros(3) zeros(3) eye(3); zeros(3) Sigma Omega];     
    
    %Controlability matrix 
    C = ctrb(A,B); 
    controlable(i) = (rank(C) == size(A,1));
end

%% GNC: LQR/SDRE control law %%
%Cost function matrices 
m = n+3;                                                    %Augmented state vector to include the integrator dynamics
Q = diag([1e-4, 1e-4, 1e-4, 1e2, 1e2, 1e2, 1e3, 1e3, 1e3]); %Cost weight on the state error
R = eye(n/2);                                               %Cost weight on the spent energy

%Preallocation 
Sc = zeros(length(tspan), 2*n);                             %Relative orbit trajectory
Sc(1,:) = Sn(1,:);                                          %Initial relative state
u = zeros(3,length(tspan));                                 %Control law
e = zeros(1,length(tspan));                                 %Error to rendezvous 

%Hamiltonian matrix of the SDRE 
H = [zeros(m,m) -B*R^(-1)*B.'; -Q zeros(m,m)];

%Initial value 
integrator = zeros(3,1);

%Compute the trajectory
for i = 1:length(tspan)
    %State coefficients 
    r_t = S(i,1:3).';                       %Synodic position of the target
        
    %Relative position between the primaries and the target 
    Ur1 = r_t-R(:,1);                       %Position of the target with respect to the first primary
    ur1 = Ur1/norm(Ur1);                    %Unit vector of the relative position of the target with respect to the primary
    Ur2 = r_t-R(:,2);                       %Position of the target with respect to the first primary
    ur2 = Ur2/norm(Ur2);                    %Unit vector of the relative position of the target with respect to the primary
        
    %Evaluate the linear model 
    Sigma = -((mup(1)/norm(Ur1)^3)+(mup(2)/norm(Ur2))^3)*eye(3)+3*((mup(1)/norm(Ur1)^3)*(ur1*ur1.')+(mup(2)/norm(Ur2)^3)*(ur2*ur2.'));
    A = [zeros(3) eye(3) zeros(3); zeros(3) zeros(3) eye(3); zeros(3) Sigma Omega];     
    
    %Compute the LQR matrix 
    H(1:9,1:9) = A;                                                              %Complete the Hamiltonian matrix
    H(end-8:end, end-8:end) = -A.';                                              %Complete the Hamiltonian matrix
    [U,L] = schur(H);                                                            %Schur decomposition of the Hamiltonian matrix
    U = U.';                                                                     %Format transformation
    P = U(1:9,10:18)*(U(1:9,1:9)^(-1)).';                                         %Solution of the Ricatti equation
    K = -R^(-1)*B.'*P;                                                           %LQR matrix
        
    %Compute the feedback control law
    u(:,i) = K*[integrator; shiftdim(Sc(i,7:12))];                                            
    
    %Re-integrate trajectory
    [~, s] = ode113(@(t,s)nlr_model(mu, true, false, 'Encke C', t, s, u(:,i)), [0 dt], S(i,:), options);
    
    %Update initial conditions
    Sc(i+1,:) = s(end,:);
    
    %Update integrator
    fintegrator = @(t,s)(shiftdim(Sc(i,7:9)));
    [~,integrator] = ode45(@(t,s)fintegrator(t,s), [0 dt], integrator, options);
    integrator = integrator(end,:).';
    
    %Error in time 
    e(i) = norm(s(end,7:12));
end

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
plot3(Sc(:,7), Sc(:,8), Sc(:,9)); 
xlabel('Synodic x coordinate');
ylabel('Synodic y coordinate');
zlabel('Synodic z coordinate');
grid on;
title('Relative motion in the configuration space');

%Configuration space error 
figure(3)
plot(tspan, e); 
xlabel('Nondimensional epoch');
ylabel('Absolute error');
grid on;
title('Absolute error in the configuration space (L2 norm)');

%Rendezvous animation 
if (false)
    figure(4) 
    view(3) 
    grid on;
    hold on
    plot3(Sc(1:index,1), Sc(1:index,2), Sc(1:index,3), 'k-.'); 
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    title('Rendezvous simulation');
    for i = 1:size(Sc,1)
        T = scatter3(Sc(i,1), Sc(i,2), Sc(i,3), 30, 'b'); 
        V = scatter3(Sc(i,1)+Sc(i,7), Sc(i,2)+Sc(i,8), Sc(i,3)+Sc(i,9), 30, 'r');

        drawnow;
        delete(T); 
        delete(V);
    end
hold off
end

