%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 10/07/21 % 

%% GNC 12: Complete rendezvous mission example 6 %% 
% This script provides an interface to test the general control scheme for a rendezvous, docking and undocking mission. 

% This mission implies the generation of user-defined artificial periodic orbits by means of low-thrust control.

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

%CR3BP constants 
mu = 3.003e-6;                      %Sun-Earth reduced gravitational parameter
Lem = 149597870700;                 %Mean distance from the Sun to the Earth

%% Initial conditions and halo orbit computation %%
%Initial conditions
L = libration_points(mu);                                   %System libration points
Az = 9000e3;                                                %Orbit amplitude out of the synodic plane. Play with it! 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 1;                                                      %Number of periods to compute
param = [1 Az Ln gamma m];                                  %Halo orbit parameters (-1 being for southern halo)

center = [0; 0; 0];                                         %Center of the artificial halo orbit  
K = 4;                                                      %Time of flight in orbital periods of the target halo
high_thrust = true;                                         %Do not use a robust control law

%Correction parameters 
dt = 1e-3;                                                  %Time step to integrate converged trajectories
maxIter = 20;                                               %Maximum allowed iterations in the differential correction schemes
tol = 1e-10;                                                %Differential correction tolerance 
Bif_tol = 1e-2;                                             %Bifucartion tolerance on the stability index
num = 1;                                                    %Number of orbits to continuate
direction = -1;                                             %Direction to continuate (to the Earth)
   
%% Functions
%Compute the halo
[halo_seed, haloT] = object_seed(mu, param, 'Halo');        %Generate a halo orbit seed

%Halo orbit 
[target_orbit, state] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Generate the initial periodic orbit
s0 = target_orbit.Trajectory(1,1:n).';
tspan = 0:dt:target_orbit.Period;
[~, S] = ode113(@(t,s)cr3bp_equations(mu, true, false, t, s), tspan, s0, options);
Sn = S; 

%Generate the target artificial periodic orbit
Sg = S;                                          %Inicialization                              
Sg(:,1:3) = Sg(:,1:3)-L(1:3,Ln).'+center.';      %Relative trajectory to the artificial Lagrange point
Sgr = Sg;                                        %Inicialization  
Sgr(:,1:3) = Sgr(:,1:3)-L(1:3,Ln).';             %Relative trajectory to the target

%% Guidance (Lissajous trajectory around the target)
%Phase definition 
tf = K*target_orbit.Period; 
tspan = 0:dt:tf; 
Sgr = repmat(Sgr(1:end-1,:),K,1);
Sgr(end+1,:) = Sgr(1,:);

%Compute the trajectory as a Chebyshev analytical expression
order = 100; 
[Cp, Cv, Cg] = CTR_guidance(order, tspan, Sgr);

%Reconstructed guidance trajectory
T = zeros(order, length(tspan));                                    %Preallocation of the polynomial basis
u = (2*tspan-(tspan(end)+tspan(1)))/(tspan(end)-tspan(1));          %Normalized time domain

for i = 1:length(tspan)
    T(:,i) = chebyshev('first', order, u(i));
end

%Error in the regression
p = Cp*T;                   %Position regression
v = Cv*T;                   %Velocity regression

%% Artificial Halo Orbit generation
if (high_thrust)
    %Artificial Halo Orbit with high thrust
    %GNC algorithms definition 
    GNC.Algorithms.Guidance = 'CTR';            %Guidance algorithm
    GNC.Algorithms.Navigation = '';             %Navigation algorithm
    GNC.Algorithms.Control = 'SMC';             %Control algorithm
    GNC.Algorithms.Solver = 'Encke';            %Dynamics vector field to be solved
    GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
    GNC.Control.Dimension = 3;                  %Dimension of the control law
    GNC.System.mu = mu;                         %System reduced gravitational parameter
    
    %Controller parameters
    %GNC.Control.SMC.Parameters = [1 SMC_optimization(mu, 'L2', St2(end,1:12), tf(2))]; 
    GNC.Control.SMC.Parameters = [1 1.000000000000000 0.432562054680836 0.070603623964497 0.099843662546135]; 
    
    %Guidance parameters 
    GNC.Guidance.CTR.Order = order;                     %Order of the approximation
    GNC.Guidance.CTR.TOF = tf;                          %Time of flight
    GNC.Guidance.CTR.PositionCoefficients = Cp;     	%Coefficients of the Chebyshev approximation
    GNC.Guidance.CTR.VelocityCoefficients = Cv;         %Coefficients of the Chebyshev approximation
    GNC.Guidance.CTR.AccelerationCoefficients = Cg;     %Coefficients of the Chebyshev approximation
     
    %Re-integrate trajectory
    s0 = [L(1:3,Ln).' 0 0 0 S(1,1:6)-[L(1:3,Ln).' 0 0 0]];
    tic
    [~, St] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
    toc
    
    %Control effort
    [~, ~, u] = GNC_handler(GNC, St(:,1:6), St(:,7:12), tspan); 
    effort = control_effort(tspan, u, false);  

else
    %Artificial Halo Orbit with low thrust
    %GNC algorithms definition 
    GNC.Algorithms.Guidance = '';               %Guidance algorithm
    GNC.Algorithms.Navigation = '';             %Navigation algorithm
    GNC.Algorithms.Control = 'TAHO';            %Control algorithm
    GNC.Algorithms.Solver = 'Encke';            %Dynamics vector field to be solved
    GNC.Guidance.Dimension = 9;                 %Dimension of the guidance law
    GNC.Control.Dimension = 3;                  %Dimension of the control law
    GNC.System.mu = mu;                         %System reduced gravitational parameter
    
    %Controller parameters
    GNC.Control.TAHO.center = center;
     
    %Re-integrate trajectory
    s0 = [L(1:3,Ln).' 0 0 0 Sgr(1,:)];
    tic
    [~, St] = ode113(@(t,s)nlr_model(mu, true, false, false, 'Encke', t, s, GNC), tspan, s0, options);
    toc
    
    %Control effort
    [~, ~, u] = GNC_handler(GNC, St(:,1:6), St(:,7:12), tspan);  
end
%% Plotting
figure(1)
view(3) 
hold on
c = plot3(Sn(:,1), Sn(:,2), Sn(:,3), 'r', 'Linewidth', 0.1); 
r = plot3(St(:,7)+St(:,1), St(:,8)+St(:,2), St(:,9)+St(:,3), 'b', 'Linewidth', 0.1); 
g = plot3(Sg(:,1), Sg(:,2), Sg(:,3), 'k', 'Linewidth', 0.1);
scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
hold off
text(L(1,Ln)+1e-3, L(2,Ln), 0, '$L_1$');
xlabel('Synodic $x$ coordinate');
ylabel('Synodic $y$ coordinate');
zlabel('Synodic $z$ coordinate');
grid on;
legend('Initial orbit', 'Rendezvous arc', 'Target orbit', 'Location', 'northeast');
title('Converged rendezvous trajectory in the absolute configuration space');

%Configuration space evolution
figure(2)
subplot(1,2,1)
hold on
plot(tspan(1:size(St,1)), St(:,7)); 
plot(tspan(1:size(St,1)), St(:,8)); 
plot(tspan(1:size(St,1)), St(:,9)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative configuration coordinate');
grid on;
legend('$x$', '$y$', '$z$');
title('Relative position evolution');

subplot(1,2,2)
hold on
plot(tspan(1:size(St,1)), St(:,10)); 
plot(tspan(1:size(St,1)), St(:,11)); 
plot(tspan(1:size(St,1)), St(:,12)); 
hold off
xlabel('Nondimensional epoch');
ylabel('Relative velocity coordinate');
grid on;
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$');
title('Relative velocity evolution');

if (false)
    dh = 50; 
    steps = fix(size(St,1)/dh);
    M = cell(1,steps);
    h = figure;
    filename = 'nhro.gif';
    view([37 20])
    hold on
    plot3(St(:,1), St(:,2), St(:,3), '.-b', 'Linewidth', 0.1);
    xlabel('Synodic x coordinate');
    ylabel('Synodic y coordinate');
    zlabel('Synodic z coordinate');
    scatter3(L(1,Ln), L(2,Ln), 0, 'k', 'filled');
    scatter3(1-mu, 0, 0, 'k', 'filled');
    text(L(1,Ln)-2e-2, L(2,Ln), 0, '$L_2$');
    text(1-mu+1e-3, 0, 5e-3, '$M_2$');
    grid on;
    title('Rendezvous simulation');
    
    for i = 1:dh:size(St,1)
        T = scatter3(St(i,1), St(i,2), St(i,3), 20, 'b', 'filled'); 
        V = scatter3(St(i,1)+St(i,7), St(i,2)+St(i,8), St(i,3)+St(i,9), 20, 'r', 'filled');
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