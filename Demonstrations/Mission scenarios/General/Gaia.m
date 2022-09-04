%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 16/05/22
% File: Gaia.m 
% Issue: 0 
% Validated: 

%% Gaia %%
% This scripts provides a test interface for the rest of the library
% functions. 

% Gaia is concerned to the generation of a gif movie of Gaia's spacecraft transfer and nominal orbit at a SEL2 Lissajous orbit.

% Credit to Grebow, 2006, for his initial seeds!

%% Test values and constants
%Set graphical environment 
set_graphics(); 

%Initial conditions
mu = 3.003e-6;                                              %Reduced gravitational parameter of the system (Sun-Earth)
Lem = 149597870700;                                         %Mean distance from the Sun to the Earth
L = libration_points(mu);                                   %System libration points
Az = 50e6;                                                  %Orbit amplitude out of the synodic plane. Play with it!
Ax = 50e6;                                                  %Orbit amplitude in the synodic plane. Play with it! 
Az = dimensionalizer(Lem, 1, 1, Az, 'Position', 0);         %Normalize distances for the E-M system
Ax = dimensionalizer(Lem, 1, 1, Ax, 'Position', 0);         %Normalize distances for the E-M system
Ln = 2;                                                     %Orbits around Li. Play with it! (L1 or L2)
gamma = L(end,Ln);                                          %Li distance to the second primary
m = 20;                                                      %Number of periods to compute
param_lyap = [Ax Az 0 0 Ln gamma m];                        %Lyapunov orbit parameters

%Correction parameters 
maxIter = 50;      %Maximum allowed iterations in the differential correction schemes
tol = 1e-10;       %Tolerance 

%Integration tolerances
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);     

%% Functions
%Compute the Lissajous seed
[halo_seed, haloT] = object_seed(mu, param_lyap, 'Lyapunov');   %Generate a Lissajous orbit

[halo_orbit, state(2)] = differential_correction('Plane Symmetric', mu, halo_seed, maxIter, tol);

%Parking orbit definition 
R = [1-mu; 0; 0];                       %Primary location in the configuration space
branch = 'L';                           %Manifold branch to globalize
map = 'Secondary primary';              %Poincaré map to use
event = @(t,s)sp_crossing(t,s,mu);      %Integration event

hd = dimensionalizer(Lem, 1, 1, 2000e3, 'Position', 0);                  %Parking orbit altitude

%Integrate the stable manifold backwards and check if it intersects the whereabouts of the parking orbit
manifold = 'S';                                                                           %Integrate the stable manifold
tspan = 0:1e-3:5*halo_orbit.Period;                                                       %Original integration time
rho = 50;                                                                                 %Density of fibres to analyze
S = invariant_manifold(mu, Ln, manifold, branch, halo_orbit.Trajectory, rho, tspan, map); %Initial trajectories

%Relative distance to the primary of interest
distance = zeros(rho,1);    
for i = 1:rho
    %Distance to the orbital altitude
    distance(i) = norm(shiftdim(S.Trajectory(i,S.ArcLength(i),1:3))-R)-hd;  
end

[~, index] = sort(distance);                            %Select the closest manifold to the parking orbit
s0 = shiftdim(S.Trajectory(index(1),1,:));              %Initial conditions to correct
Phi = eye(n);                                           %Initial STM 
Phi = reshape(Phi, [n^2 1]);                            %Initial STM 
sHalo = seed(S.Index(index(1)),1:n).';                  %Halo insertion point
s0 = [s0(1:3); s0(4:6); Phi];                           %Initial conditions
TOF = S.TOF(index(1));                                  %Time of flight
tspan = TOF:-dt:0;                                                                  %Integration time
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, 'Events', event);             %Integration tolerances  
[~, St0] = ode113(@(t,s)cr3bp_equations(mu, 1, true, t, s), tspan, s0, options);    %Natural trajectory

%Complete trajectory 
Stotal = [St0; halo_seed];

%% Plotting
if (true)
    dh = 100;
    W = figure(1);
    filename = 'gaia.gif';
    view([40 20]) 
    hold on
    H = plot3(halo_seed(:,1), halo_seed(:,2), halo_seed(:,3), 'b');
    H.Color(4) = 0.3;
    for i = 1:size(S.Trajectory,1)
        ManifoldAux = shiftdim(S.Trajectory(i,:,:));
        plot3(ManifoldAux(1:S.ArcLength(i),1), ManifoldAux(1:S.ArcLength(i),2), ManifoldAux(1:S.ArcLength(i),3), 'g');
    end
    plot3(St0(1:2000,1), St0(1:2000,2), St0(1:2000,3), 'm');
    scatter3(L(1,1), L(2,1), 0, 'k', 'filled');
    scatter3(L(1,2), L(2,2), 0, 'k', 'filled');
    scatter3(1-mu, 0, 0, 'k', 'filled');
    text(L(1,1)+1e-3, L(2,1)+1e-3, 5e-3, '$L_1$');
    text(L(1,2)+1e-3, L(2,2), 5e-3, '$L_2$');
    text(1-mu+1e-3, 0, 5e-3, '$M_2$');
    xlabel('Synodic $x$ coordinate');
    ylabel('Synodic $y$ coordinate');
    zlabel('Synodic $z$ coordinate');
    grid on;
    legend('Gaia Lissajous orbit', 'Transfer trajectory', 'Stable manifold', 'Location', 'northeast', 'AutoUpdate', 'off');
    title('Rendezvous simulation');

    for i = 1:dh:size(Stotal,1)
        T = scatter3(Stotal(i,1), Stotal(i,2), Stotal(i,3), 30, 'b', 'filled');

        drawnow;
        frame = getframe(W);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); 
        if (i == 1) 
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1e-3); 
        else 
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1e-3); 
        end 
        delete(T); 
    end
    hold off
end