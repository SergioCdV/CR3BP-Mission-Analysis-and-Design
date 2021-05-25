%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 01/04/21 % 

%% GNC 3: LQR/SDRE control law %% 
% This script provides an interface to test LQR rendezvous strategies for
% rendezvous missions.

% The relative motion of two spacecraft in the same halo orbit (closing and RVD phase) around L1 in the
% Earth-Moon system is analyzed.

% The classical LQR and the SDRE, with their discrete versions, in three
% different linear models.

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
tf = 2*pi;                          %Rendezvous time
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

%% Controlability analysis
model = 'Target';
controlable = controlability_analysis(model, mu, S, index, Ln, gamma);

%% GNC algorithms definition 
GNC.Algorithms.Guidance = '';                   %Guidance algorithm
GNC.Algorithms.Navigation = '';                 %Navigation algorithm
GNC.Algorithms.Control = 'SDRE';                 %Control algorithm

GNC.Guidance.Dimension = 9;                     %Dimension of the guidance law
GNC.Control.Dimension = 3;                      %Dimension of the control law

GNC.System.mu = mu;                             %Systems's reduced gravitational parameter
GNC.System.Libration = [Ln gamma];              %Libration point ID

GNC.Control.LQR.Model = model;                  %LQR model
GNC.Control.SDRE.Model = model;                 %SDRE model
GNC.Control.LQR.Q = eye(9);                     %Penalty on the state error
GNC.Control.LQR.M = eye(3);                     %Penalty on the control effort
GNC.Control.LQR.Reference = Sn(index,1:3);      %Penalty on the control effort
GNC.Control.SDRE.Q = eye(9);                    %Penalty on the state error
GNC.Control.SDRE.M = eye(3);                    %Penalty on the control effort

%% GNC: SDRE/LQR control law
%Initial conditions 
int = zeros(1,3);                               %Integral of the relative position
slqr0 = [Sn(1,:) int];                          %Initial conditions

%Compute the trajectory
[~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, slqr0, options);

%Error in time 
[e, merit] = figures_merit(tspan, Sc);

%Control law
[~, ~, u] = GNC_handler(GNC, Sc(:,1:6), Sc(:,7:end), NaN);

%Control integrals
energy = control_effort(tspan, u);

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

%Configuration space evolution
figure(3)
subplot(1,2,1)
hold on
plot(tspan, Sc(:,7)); 
plot(tspan, Sc(:,8)); 
plot(tspan, Sc(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('x coordinate', 'y coordinate', 'z coordinate');
title('Relative position evolution');
subplot(1,2,2)
hold on
plot(tspan, Sc(:,10)); 
plot(tspan, Sc(:,11)); 
plot(tspan, Sc(:,12)); 
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
ylabel('Absolute error (log)');
grid on;
title('Absolute error in the configuration space (L2 norm)');

%Rendezvous animation 
if (false)
    figure(5) 
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

%% Auxiliary functions
%Controlability analysis function 
function [controlable] = controlability_analysis(model, mu, Sn, index, Ln, gamma)
    %Approximation 
    n = 6;                              %Dimension of the state vector
    order = 2;                          %Order of the approximation 

    %Preallocation 
    controlable = zeros(index,1);       %Controllability boolean

    %Model coefficients 
    mup(1) = 1-mu;                      %Reduced gravitational parameter of the first primary 
    mup(2) = mu;                        %Reduced gravitational parameter of the second primary 
    R(:,1) = [-mu; 0; 0];               %Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];              %Synodic position of the second primary

    %Linear model matrices
    B = [zeros(n/2); zeros(n/2); eye(n/2)];         %Linear model input matrix 
    Omega = [0 2 0; -2 0 0; 0 0 0];                 %Coriolis dyadic

    for i = 1:index
        %State coefficients 
        r_t = Sn(i,1:3).';                          %Synodic position of the target
        
        %Select linear model 
        switch (model)
            case 'Steady libration point'
                cn = legendre_coefficients(mu, Ln, gamma, order);     %Compute the relative Legendre coefficient c2 
                c2 = cn(2); 
                Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              %Linear model state matrix
            case 'Unsteady libration point' 
                cn = relegendre_coefficients(mu, r_t.', order);       %Compute the relative Legendre coefficient c2 
                c2 = cn(2); 
                Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              %Linear model state matrix
            case 'Target' 
                %Relative position between the primaries and the target 
                Ur1 = r_t-R(:,1);               %Position of the target with respect to the first primary
                ur1 = Ur1/norm(Ur1);            %Unit vector of the relative position of the target with respect to the primary
                Ur2 = r_t-R(:,2);               %Position of the target with respect to the first primary
                ur2 = Ur2/norm(Ur2);            %Unit vector of the relative position of the target with respect to the primary

                %Evaluate the linear model 
                Sigma = -((mup(1)/norm(Ur1)^3)+(mup(2)/norm(Ur2))^3)*eye(3)+3*((mup(1)/norm(Ur1)^3)*(ur1*ur1.')+(mup(2)/norm(Ur2)^3)*(ur2*ur2.'));
            otherwise 
                error('No valid linear model was selected'); 
        end

        %Linear state model
        A = [zeros(3) eye(3) zeros(3); zeros(3) zeros(3) eye(3); zeros(3) Sigma Omega];     

        %Controlability matrix 
        C = ctrb(A,B); 
        controlable(i) = (rank(C) == size(A,1));
    end
end

%Implement a LQR control law 
function [Sc, e] = lqrm(model, options, mu, Sn, tspan, Ln, gamma)
    %Approximation 
    n = 6;                              %Dimension of the state vector
    order = 2;                          %Order of the approximation 

    %Model coefficients 
    mup(1) = 1-mu;                      %Reduced gravitational parameter of the first primary 
    mup(2) = mu;                        %Reduced gravitational parameter of the second primary 
    R(:,1) = [-mu; 0; 0];               %Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];              %Synodic position of the second primary

    %Linear model matrices
    B = [zeros(n/2); zeros(n/2); eye(n/2)];         %Linear model input matrix 
    Omega = [0 2 0; -2 0 0; 0 0 0];                 %Coriolis dyadic
    
    %Cost function matrices 
    Q = diag(ones(1,n+3));                          %Cost weight on the state error
    M = eye(n/2);                                   %Cost weight on the spent energy

    %Preallocation 
    e = zeros(1,length(tspan));                     %Error to rendezvous 

    %Initial conditions 
    int = zeros(1,3);
    slqr0 = [Sn(1,:) int];

    %State coefficients 
    r_t = Sn(end,1:3).';                              %Synodic position of the target

    %Select linear model 
    switch (model)
        case 'Steady libration point'
            cn = legendre_coefficients(mu, Ln, gamma, order);     %Compute the relative Legendre coefficient c2 
            c2 = cn(2); 
            Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              %Linear model state matrix
        case 'Unsteady libration point' 
            rc = relegendre_coefficients(mu, r_t.', order);       %Compute the relative Legendre coefficient c2 
            c2 = rc(2); 
            Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              %Linear model state matrix
        case 'Target' 
            %Relative position between the primaries and the target 
            Ur1 = r_t-R(:,1);               %Position of the target with respect to the first primary
            ur1 = Ur1/norm(Ur1);            %Unit vector of the relative position of the target with respect to the primary
            Ur2 = r_t-R(:,2);               %Position of the target with respect to the first primary
            ur2 = Ur2/norm(Ur2);            %Unit vector of the relative position of the target with respect to the primary
            %Evaluate the linear model 
            Sigma = -((mup(1)/norm(Ur1)^3)+(mup(2)/norm(Ur2))^3)*eye(3)+3*((mup(1)/norm(Ur1)^3)*(ur1*ur1.')+(mup(2)/norm(Ur2)^3)*(ur2*ur2.'));
        otherwise 
            error('No valid linear model was selected'); 
    end

    %Linear state model
    A = [zeros(3) eye(3) zeros(3); zeros(3) zeros(3) eye(3); zeros(3) Sigma Omega];  

    %Compute the feedback control law
    [K,~,~] = lqr(A,B,Q,M);

    %Compute the trajectory
    [~, Sc] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke LQR', t, s, K), tspan, slqr0, options);

    %Error in time 
    for i = 1:length(tspan)
        e(i) = norm(Sc(i,7:12));
    end
end
