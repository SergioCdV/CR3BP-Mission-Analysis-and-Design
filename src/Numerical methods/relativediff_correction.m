%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 03/01/21
% File: differential_correction.m 
% Issue: 0 
% Validated: 

%% Relative Differential correction %%
% This function contains the algorithm to compute several correction
% algorithms depending on some user inputs.

% Inputs: - string algorithm, selecting the differential correction scheme to use. 
%         - double mu, the reduced gravitational parameter of the system.
%         - object seed, which vary depending on the algorithm in use. 
%         - int maxIter, number of maximum allowed corrections.
%         - double tol, tolerance to stop the correction. 

% Outputs: - vector xf, converged/last corrected trajectory.
%          - boolean state, outputing the convergence of the algorithm.

% Methods: multiple shooting is employed for general periodic orbits if a full trajectory seed is available.
%          Single shooting with symmetric constraints is employed for computing relative Lyapunov and Halo orbits, 
%          although convergence may be worse. Torus invariant curves are corrected via the method by Scheers and Haapala.

% New versions: torus correction. Use monodromy properties to correct error.

function [xf, state] = relativediff_correction(algorithm, mu, seed, maxIter, tol, varargin)            
    %Implement the selected scheme 
    switch (algorithm)
        case 'Axis Symmetric'
            [xf, state] = Sym_Axis_scheme(mu, seed, maxIter, tol);
        case 'Plane Symmetric'
            [xf, state] = Sym_Plane_scheme(mu, seed, maxIter, tol);
        case 'Double Symmetric'
            [xf, state] = Sym_Double_scheme(mu, seed, maxIter, tol);
        case 'Double Plane Symmetric' 
            [xf, state] = Sym_DoublePlane_scheme(mu, seed, maxIter, tol);
        case 'Planar'
            [xf, state] = Sym_Planar_scheme(mu, seed, maxIter, tol);
        otherwise
            error('No valid option was selected');
    end
end

%% Auxiliary functions (individual schemes)
%Compute periodic orbits using the X axis symmetry
function [xf, state] = Sym_Axis_scheme(mu, seed, maxIter, tol) 
    %Constants 
    m = 12;      %Phase space dimension 
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == m) || (size(seed,1) == m)
        %Restrict the seed to the initial conditions
        if (size(seed,2) == m)
            if (size(seed,1) ~= 1)
                seed = shiftdim(seed(1,:));
            else
                seed = seed.';
            end
        elseif (size(seed,1) == m) 
            if (size(seed,2) ~= 1)
                seed = shiftdim(seed(:,1));
            end
        end
    else
        error('No valid initial conditions were input');
    end
    
    %Ensure required initial conditions
    seed(m/2+2) = 0;    %Null Y coordinate 
    seed(m/2+3) = 0;    %Null Z coordinate
    seed(m/2+4) = 0;    %Null Vx 
    
    %Augment initial conditions with the initial STM 
    Phi = eye(m/2); 
    Phi = reshape(Phi, [(m/2)^2 1]); 
    seed(m+1:m+(m/2)^2) = Phi;
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, ... 
                     'Events', @(t,s)x_crossing(t,s));          %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    T = 2*pi;                                                   %Initial orbit period
    tspan = 0:dt:T;                                             %Integration time span
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;                %Convergence flag
    iter = 1;                   %Initial iteration
    
    %Preallocation 
    ds0 = zeros(3,maxIter);     %Vector containing the initial conditions correction
        
    %Main computation 
    while (GoOn) && (iter < maxIter)
        %Proceed with the integration
        [~, S] = ode113(@(t,s)nlr_model(mu, direction, false, flagVar, 'Encke', t, s), tspan, seed, options);
        
        %Compute error
        e = [S(end,m/2+3); S(end,m/2+4)];   %Z and Vx at the crossing must be 0
        
        %Compute the correction 
        F = nlr_model(mu, direction, false, flagVar, 'Encke', 0, S(end,:).');         %Vector field at T/2
        Phi = reshape(S(end,m+1:end), [m/2 m/2]);                                     %Build the monodromy matrix at T/2
        A = [Phi(3:4,1) Phi(3:4,5) Phi(3:4,6)] ...
            -(1/S(end,m/2+5)*[S(end,m/2+6); F(4)]*[Phi(2,1) Phi(2,5) Phi(2,6)]);      %Constraint matrix
        ds0(:,iter) = pinv(A)*e;                                                      %Compute the variation
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            seed(m/2+1) = seed(m/2+1)-ds0(1,iter);                  %Update initial x coordinate
            seed(m/2+5:m/2+6) = seed(m/2+5:m/2+6)-ds0(2:3,iter);    %Update initial conditions
            iter = iter+1;                                          %Update iteration
        end       
    end
    
    %Integrate the whole trayectory 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Disable crossing event
    tspan = 0:dt:(2*dt*size(S,1));                              %Integrate the orbit for a whole period
    
    [t, S] = ode113(@(t,s)nlr_model(mu, direction, false, flagVar, 'Encke', t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;          %Trajectory
    xf.Period = t(end);         %Orbit period
    
    %Ouput differential correction scheme convergence results
    state.State = ~GoOn;        %Convergence boolean
    state.Iterations = iter;    %Final iteration
    state.Error = norm(e);      %Final error L2 norm
end

%Compute periodic orbits using the XZ symmetry
function [xf, state] = Sym_Plane_scheme(mu, seed, maxIter, tol)
    %Constants 
    m = 12;                  %Phase space dimension 

    %Sanity check on initial conditions dimension
    if (size(seed,2) == m) || (size(seed,1) == m)
        %Restrict the seed to the initial conditions
        if (size(seed,2) == m)
            if (size(seed,1) ~= 1)
                seed = shiftdim(seed(1,:));
            else
                seed = seed.';
            end
        elseif (size(seed,1) == m) 
            if (size(seed,2) ~= 1)
                seed = shiftdim(seed(:,1));
            end
        end
    else
        error('No valid initial conditions were input');
    end
    
    %Ensure required initial conditions
    seed(m/2+2) = 0;    %Null Y coordinate 
    seed(m/2+4) = 0;    %Null Vx 
    seed(m/2+6) = 0;    %Null Vz
    
    %Augment initial conditions with the initial STM 
    Phi = eye(m/2); 
    Phi = reshape(Phi, [(m/2)^2 1]); 
    seed(m+1:m+(m/2)^2) = Phi;
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, ... 
                     'Events', @(t,s)x_crossing(t,s));          %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    T = 2*pi;                                                   %Initial orbit period
    tspan = 0:dt:T;                                             %Integration time span
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;                %Convergence flag 
    iter = 1;                   %Initial iteration
    
    %Preallocation 
    ds0 = zeros(2,maxIter);     %Vector containing the initial conditions correction

    %Main computation 
    while (GoOn) && (iter < maxIter)
        %Proceed with the integration
        [~, S] = ode113(@(t,s)nlr_model(mu, direction, false, flagVar, 'Encke', t, s), tspan, seed, options);
        
        %Compute error
        e = [S(end,m/2+4); S(end,m/2+6)];   %Vx and Vz at the crossing must be 0
        
        %Compute the correction 
        F = nlr_model(mu, direction, false, flagVar, 'Encke', 0, S(end,:).');         %Vector field at T/2
        Phi = reshape(S(end,m+1:end), [m/2 m/2]);                                     %Build the monodromy matrix at T/2
        A = [Phi(4,1) Phi(4,5); Phi(6,1) Phi(6,5)] ...
            -(1/S(end,m/2+5)*[F(4); F(6)]*[Phi(2,1) Phi(2,5)]);                       %Constraint matrix
        
        %Update the initial conditions
        ds0(:,iter) = pinv(A)*e;
        seed(m/2+1) = seed(m/2+1)-ds0(1,iter);            %Update the initial x coordinate
        seed(m/2+5) = seed(m/2+5)-ds0(2,iter);            %Update the initial Vy coordinate
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            iter = iter+1;                        %Update iteration
        end       
    end
    
    %Integrate the whole trayectory 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Disable crossing event
    tspan = 0:dt:(2*dt*size(S,1));                              %Integrate the orbit for a whole orbit
    [t, S] = ode113(@(t,s)nlr_model(mu, direction, false, flagVar, 'Encke', t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;          %Trajectory
    xf.Period = t(end);         %Orbit period
    
    %Ouput differential correction scheme convergence results
    state.State = ~GoOn;        %Convergence boolean
    state.Iterations = iter;    %Final iteration
    state.Error = norm(e);      %Final error L2 norm
end

%Compute periodic orbits using the double X-XZ symmetry
function [xf, state] = Sym_Double_scheme(mu, seed, maxIter, tol)
    %Constants 
    m = 12;      %Phase space dimension 
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == m) || (size(seed,1) == m)
        %Restrict the seed to the initial conditions
        if (size(seed,2) == m)
            if (size(seed,1) ~= 1)
                seed = shiftdim(seed(1,:));
            else
                seed = seed.';
            end
        elseif (size(seed,1) == m) 
            if (size(seed,2) ~= 1)
                seed = shiftdim(seed(:,1));
            end
        end
    else
        error('No valid initial conditions were input');
    end
    
    %Ensure required initial conditions
    seed(m/2+2) = 0;    %Null Y coordinate 
    seed(m/2+3) = 0;    %Null Z coordinate
    seed(m/2+4) = 0;    %Null Vx 
    
    %Augment initial conditions with the initial STM 
    Phi = eye(m/2); 
    Phi = reshape(Phi, [(m/2)^2 1]); 
    seed(m+1:m+(m/2)^2) = Phi;
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, ... 
                     'Events', @(t,s)x_crossing(t,s));          %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    T = 2*pi;                                                   %Initial orbit period
    tspan = 0:dt:T;                                             %Integration time span
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;                %Convergence flag   
    iter = 1;                   %Initial iteration
    
    %Preallocation 
    ds0 = zeros(3,maxIter);     %Vector containing the initial conditions correction
        
    %Main computation 
    while (GoOn) && (iter < maxIter)
        %Proceed with the integration
        [~, S] = ode113(@(t,s)nlr_model(mu, direction, false, flagVar, 'Encke', t, s), tspan, seed, options);
        
        %Compute error
        e = [S(end,m/2+4); S(end,m/2+6)];   %Vx and Vz at the crossing must be 0
        
        %Compute the correction 
        F = nlr_model(mu, direction, false, flagVar, 'Encke', 0, S(end,:).');         %Vector field at T/4
        Phi = reshape(S(end,m+1:end), [m/2 m/2]);                                     %Build the monodromy matrix at T/4
        A = [Phi(4,1) Phi(4,5) Phi(4,6); Phi(6,1) Phi(6,5) Phi(6,6)] ...
            -(1/S(end,m/2+5)*[F(4); F(6)]*[Phi(2,1) Phi(2,5) Phi(2,6)]);              %Constraint matrix
        ds0(:,iter) = pinv(A)*e;                                                      %Compute the variation
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            seed(m/2+1) = seed(m/2+1)-ds0(1,iter);                  %Update initial x coordinate
            seed(m/2+5:m/2+6) = seed(m/2+5:m/2+6)-ds0(2:3,iter);    %Update initial conditions
            iter = iter+1;                                          %Update iteration
        end       
    end
    
    %Integrate the whole trayectory 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Disable crossing event
    tspan = 0:dt:(4*dt*size(S,1));                              %Integrate the orbit for a whole period
    [t, S] = ode113(@(t,s)nlr_model(mu, direction, false, flagVar, 'Encke', t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;          %Trajectory
    xf.Period = t(end);         %Orbit period
    
    %Ouput differential correction scheme convergence results
    state.State = ~GoOn;        %Convergence boolean
    state.Iterations = iter;    %Final iteration
    state.Error = norm(e);      %Final error L2 norm
end

%Compute periodic orbits using the symmetry XZ-XY planes
function [xf, state] = Sym_DoublePlane_scheme(mu, seed, maxIter, tol)
    %Constants 
    m = 12;      %Phase space dimension 
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == 12) || (size(seed,1) == 12)
        %Restrict the seed to the initial conditions
        if (size(seed,2) == 12)
            if (size(seed,1) ~= 1)
                seed = shiftdim(seed(1,:));
            else
                seed = seed.';
            end
        elseif (size(seed,1) == 12) 
            if (size(seed,2) ~= 1)
                seed = shiftdim(seed(:,1));
            end
        end
    else
        error('No valid initial conditions were input');
    end
    
    %Ensure required initial conditions
    seed(2) = 0;    %Null Y coordinate 
    seed(4) = 0;    %Null Vx 
    seed(6) = 0;    %Null Vz
    
    %Augment initial conditions with the initial STM 
    Phi = eye(m/2); 
    Phi = reshape(Phi, [(m/2)^2 1]); 
    seed(m+1:m+(m/2)^2) = Phi;
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, ... 
                     'Events', @(t,s)x_crossing(t,s));          %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    T = 2*pi;                                                   %Initial orbit period
    tspan = 0:dt:T;                                             %Integration time span
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;                %Convergence flag 
    iter = 1;                   %Initial iteration
    
    %Preallocation 
    ds0 = zeros(2,maxIter);     %Vector containing the initial conditions correction
        
    %Main computation 
    while (GoOn) && (iter < maxIter)
        %Proceed with the integration
        [~, S] = ode113(@(t,s)nlr_model(mu, direction, false, flagVar, 'Encke', t, s), tspan, seed, options);
        
        %Compute error
        e = [S(end,m/2+4); S(end,m/2+6)];   %Vx and Vz at the crossing must be 0
        
        %Compute the correction 
        F = nlr_model(mu, direction, false, flagVar, 'Encke', 0, S(end,:).');         %Vector field at T/2
        Phi = reshape(S(end,m+1:end), [m/2 m/2]);                                     %Build the monodromy matrix at T/2
        A = [Phi(4,1) Phi(4,6); Phi(6,1) Phi(6,6)] ...
            -(1/S(end,m/2+5)*[F(4); F(6)]*[Phi(2,1) Phi(2,6)]);
        ds0(:,iter) = A\e;                                                            %Compute the variation
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            seed(m/2+1) = seed(m/2+1)-ds0(1,iter);        %Update initial conditions
            seed(m/2+6) = seed(m/2+6)-ds0(2,iter);        %Update initial conditions
            iter = iter+1;                                %Update iteration
        end       
    end
    
    %Integrate the whole trayectory 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Disable crossing event
    tspan = 0:dt:(2*dt*size(S,1));                              %Integrate the orbit for a whole orbit
    [t, S] = ode113(@(t,s)nlr_model(mu, direction, false, flagVar, 'Encke', t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;          %Trajectory
    xf.Period = t(end);         %Orbit period
    
    %Ouput differential correction scheme convergence results
    state.State = ~GoOn;        %Convergence boolean
    state.Iterations = iter;    %Final iteration
    state.Error = norm(e);      %Final error L2 norm
end

%Compute planar periodic orbits -for Lyapunov orbits-
function [xf, state] = Sym_Planar_scheme(mu, seed, maxIter, tol)
    %Constants 
    m = 12;      %Phase space dimension 
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == m) || (size(seed,1) == m)
        %Restrict the seed to the initial conditions
        if (size(seed,2) == m)
            if (size(seed,1) ~= 1)
                seed = shiftdim(seed(1,:));
            else
                seed = seed.';
            end
        elseif (size(seed,1) == m) 
            if (size(seed,2) ~= 1)
                seed = shiftdim(seed(:,1));
            end
        end
    else
        disp('No valid initial conditions were input');
    end
    
    %Ensure motion on the XY synodic plane
    seed(m/2+2) = 0;    %Null Y coordinate 
    seed(m/2+3) = 0;    %Null Z coordinate
    seed(m/2+4) = 0;    %Null Vx 
    seed(m/2+6) = 0;    %Null Vz 
    
    %Augment initial conditions with the initial STM 
    Phi = eye(m/2); 
    Phi = reshape(Phi, [(m/2)^2 1]); 
    seed(m+1:m+(m/2)^2) = Phi;
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, ... 
                     'Events', @(t,s)x_crossing(t,s));          %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    T = 2*pi;                                                   %Initial orbit period
    tspan = 0:dt:T;                                             %Integration time span
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;                %Convergence flag 
    iter = 1;                   %Initial iteration
    
    %Preallocation 
    ds0 = zeros(1,maxIter);     %Vector containing the initial conditions correction
        
    %Main computation 
    while (GoOn) && (iter < maxIter)
        %Proceed with the integration
        [~, S] = ode113(@(t,s)nlr_model(mu, direction, false, flagVar, 'Encke', t, s), tspan, seed, options);
        
        %Compute error
        e = S(end,m/2+4);   %Vx at the crossing must be 0
        
        %Compute correction 
        F = nlr_model(mu, direction, false, flagVar, 'Encke', 0, S(end,:).');   %Vector field at T/2
        Phi = reshape(S(end,m+1:end), [m/2 m/2]);                               %Build the monodromy matrix at T/2
        ds0(iter) = e/(Phi(4,5)-(F(4)/S(end,m/2+5))*Phi(2,5));                  %Compute the variation
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            seed(m/2+5) = seed(m/2+5)-ds0(iter);    %Update initial conditions
            iter = iter+1;                  %Update iteration
        end       
    end
    
    %Integrate the whole trayectory 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Disable crossing event
    tspan = 0:dt:(2*dt*size(S,1));
    [t, S] = ode113(@(t,s)nlr_model(mu, direction, false, flagVar, 'Encke', t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;          %Trajectory
    xf.Period = t(end);         %Orbit period
    
    %Ouput differential correction scheme convergence results
    state.State = ~GoOn;        %Convergence boolean
    state.Iterations = iter;    %Final iteration
    state.Error = norm(e);      %Final error L2 norm
end