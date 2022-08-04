%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 25/07/22 % 

%% GNC 15: Impulsive Center Manifold Phasing %% 
% This script provides an interface to test the ICP guidance core.

% The relative motion of two spacecraft in the two different halo orbit is analysed 
% and a guidance rendezvous trajectory is designed.

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
Az = 10e6;                                                          %Orbit amplitude out of the synodic plane. 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);                 %Normalize distances for the E-M system
Ln = 1;                                                             %Orbits around L1
gamma = L(end,Ln);                                                  %Li distance to the second primary
m = 1;                                                              %Number of periods to compute

%Compute a halo seed 
halo_param = [1 Az Ln gamma m];                                     %Northern halo parameters
[halo_seed, period] = object_seed(mu, halo_param, 'Halo');          %Generate a halo orbit seed

%Correct the seed and obtain initial conditions for a halo orbit
[target_orbit, ~] = differential_correction('Planar', mu, halo_seed, maxIter, tol);
chaser_orbit = target_orbit;

%% Modelling in the synodic frame %%
index = 100;                                                        %Phasing point
r_t0 = target_orbit.Trajectory(index,1:6);                          %Initial target conditions
r_c0 = chaser_orbit.Trajectory(1,1:6);                              %Initial chaser conditions 
rho0 = r_c0-r_t0;                                                   %Initial relative conditions
s0 = [r_t0 rho0].';                                                 %Initial conditions of the target and the relative state

% Phasing 
target_orbit.Trajectory = [target_orbit.Trajectory(index:end,:); target_orbit.Trajectory(1:index,:)];

%% Generate the guidance trajectory
% Guidance trajectory parameters
k = 10;                                         % Number of phasing revolutions
dtheta = 2*pi/target_orbit.Period*(index*dt);   % Initial phase difference

%Absolute tori
% L0 = [L(1:3,1).' zeros(1,3)];
% [Str, state(1)] = differential_rtorus(mu, target_orbit.Period, [L0 r_c0], tol);

% Phasing tori 
[Str, state(1)] = ICP_guidance(mu, target_orbit.Period, dtheta, k, [r_t0 r_c0], tol);


phasing_orbit.Trajectory = [];
phasing_orbit.Period = Str.Period;
for i = 1:size(Str.Trajectory,1)
    % Integration of the quasi-periodic trajectory
    tspan = 0:dt:Str.Period;
    [~, St] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s), tspan, Str.Trajectory(i,:), options);
    S = St(:,1:6)+St(:,7:12);
    phasing_orbit.Trajectory = [phasing_orbit.Trajectory; S(1:end-1,:)];
end 

% Periodicity check 
index = floor(mod(k*Str.Period,target_orbit.Period)/target_orbit.Period*size(target_orbit.Trajectory,1));
chaser_orbit.Trajectory = repmat(chaser_orbit.Trajectory(1:end-1,:), floor(k*Str.Period/target_orbit.Period), 1);
target_orbit.Trajectory = repmat(target_orbit.Trajectory(1:end-1,:), floor(k*Str.Period/target_orbit.Period), 1); 
target_orbit.Trajectory = [target_orbit.Trajectory; target_orbit.Trajectory(1:index,:)];

%% Results 
%Plot results 
figure(1) 
view(3) 
hold on
T = plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), 'r'); 
plot3(chaser_orbit.Trajectory(1,1), chaser_orbit.Trajectory(1,2),chaser_orbit.Trajectory(1,3), '*k')
plot3(phasing_orbit.Trajectory(:,1), phasing_orbit.Trajectory(:,2), phasing_orbit.Trajectory(:,3), 'b', 'Linewidth', 0.05)
legend('Reference target orbit', 'Chaser orbit', 'Guidance orbit', 'AutoUpdate', 'off')
plot3(S(1,1), S(1,2), S(1,3), '*r');
S0 = [target_orbit.Trajectory(1,1:6) zeros(1,6) reshape(eye(6), [1 36])]; 
plot3(target_orbit.Trajectory(1,1), target_orbit.Trajectory(1,2), target_orbit.Trajectory(1,3), '*b'); 
plot3(target_orbit.Trajectory(end,1), target_orbit.Trajectory(end,2), target_orbit.Trajectory(end,3), '*b'); 
hold off
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
title('Guidance trajectory between periodic orbits');

%Configuration space evolution
figure(2)
subplot(1,2,1)
hold on
plot(tspan(1:size(St,1)), St(:,7)); 
plot(tspan(1:size(St,1)), St(:,8)); 
plot(tspan(1:size(St,1)), St(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Configuration coordinates');
grid on;
legend('$x$', '$y$', '$z$');
title('Position in time');
subplot(1,2,2)
hold on
plot(tspan(1:size(St,1)), St(:,10)); 
plot(tspan(1:size(St,1)), St(:,11)); 
plot(tspan(1:size(St,1)), St(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Velocity coordinates');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Velocity in time');

%Rendezvous animation 
if (false)
    dh = 50; 
    steps = fix(size(St,1)/dh);
    M = cell(1,steps);
    h = figure;
    filename = 'phasing_tori.gif';
    view([37 20])
    hold on
    plot3(target_orbit.Trajectory(:,1), target_orbit.Trajectory(:,2), target_orbit.Trajectory(:,3), '.-b', 'Linewidth', 0.1);
    plot3(phasing_orbit.Trajectory(:,1), phasing_orbit.Trajectory(:,2), phasing_orbit.Trajectory(:,3), '.-b', 'Linewidth', 0.1);
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
    text(L(1,Ln)-2e-2, L(2,Ln), 0, '$L_2$');
    grid on;
    title('Rendezvous simulation');
    
    for i = 1:dh:size(St,1)
        T = scatter3(target_orbit.Trajectory(i,1), target_orbit.Trajectory(i,2), target_orbit.Trajectory(i,3), 20, 'b', 'filled'); 
        V = scatter3(chaser_orbit.Trajectory(i,1), chaser_orbit.Trajectory(i,2), chaser_orbit.Trajectory(i,3), 20, 'b', 'filled'); 
        C = scatter3(phasing_orbit.Trajectory(i,1), phasing_orbit.Trajectory(i,2), phasing_orbit.Trajectory(i,3), 20, 'r', 'filled'); 
        drawnow;
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); 
        if (i == 1) 
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1e-3); 
        else 
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1e-3); 
        end 
        delete(T); 
        delete(V);
    end
    hold off
end