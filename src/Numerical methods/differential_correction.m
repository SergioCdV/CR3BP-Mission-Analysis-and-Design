%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 03/01/21
% File: differential_correction.m 
% Issue: 0 
% Validated: 

%% Differential correction %%
% This function contains the algorithm to compute several correction
% algorithms depending on some user inputs.

% Inputs: - string algorithm, selecting the differential correction scheme to use. 
%         - double mu, the reduced gravitational parameter of the system.
%         - object seed, which vary depending on the algorithm in use. 
%         - int n, number of maximum allowed corrections.
%         - double tol, tolerance to stop the correction. 

% Outputs: - vector xf, converged/last corrected trajectory.
%          - boolean state, outputing the convergence of the algorithm.

% Methods: multiple shooting is employed for general periodic orbits if a full trajectory seed is available.
%          Single shooting with symmetric constraints is employed for computing Lyapunov and Halo orbits, 
%          although convergence may be worse. Torus invariant curves are corrected via the method by Scheers and Haapala.

% New versions: torus correction. Use monodromy properties to correct
% error. Review all.

function [xf, state] = differential_correction(algorithm, mu, seed, n, tol, varargin)            
    %Implement the selected scheme 
    switch (algorithm)
        case 'Axis Symmetric'
            [xf, state] = SymAxis_scheme(mu, seed, n, tol);
        case 'Plane Symmetric'
            [xf, state] = SymPlane_scheme(mu, seed, n, tol);
        case 'Double Symmetric'
            [xf, state] = SymDouble_scheme(mu, seed, n, tol);
        case 'Double Plane Symmetric' 
            [xf, state] = SymDoublePlane_scheme(mu, seed, n, tol);
        case 'Planar'
            [xf, state] = SymPlanar_scheme(mu, seed, n, tol);
        case 'Periodic Multiple Shooting'
            [xf, state]= MSPeriodic_scheme(mu, seed, n, tol, varargin);
        case 'Jacobi Constant Multiple Shooting'
            [xf, state] = MSJacobi_scheme(mu, seed, n, tol, varargin);
        case 'Periodic PSC Multiple Shooting'
            [xf, state] = PSPeriodic_scheme(mu, seed, n, tol, varargin);
        otherwise
            disp('No valid options was selected.');
            xf = [];
            state = false;
    end
end

%% Auxiliary functions (individual schemes)
%Compute periodic orbits using the X axis symmetry
function [xf, state] = SymAxis_scheme(mu, seed, n, tol) 
    %Constants 
    m = 6;      %Phase space dimension 
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == 6) || (size(seed,1) == 6)
        %Restrict the seed to the initial conditions
        if (size(seed,2) == 6)
            if (size(seed,1) ~= 1)
                seed = shiftdim(seed(1,:));
            else
                seed = seed.';
            end
        elseif (size(seed,1) == 6) 
            if (size(seed,2) ~= 1)
                seed = shiftdim(seed(:,1));
            end
        end
    else
        disp('No valid initial conditions.');
        xf = []; 
        state = false; 
        return;
    end
    
    %Ensure required initial conditions
    seed(2) = 0;    %Null Y coordinate 
    seed(3) = 0;    %Null Z coordinate
    seed(4) = 0;    %Null Vx 
    
    %Augment initial conditions with the initial STM 
    Phi = eye(m); 
    Phi = reshape(Phi, [m^2 1]); 
    seed(m+1:m+m^2) = Phi;
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, ... 
                     'Events', @(t,s)x_crossing(t,s));          %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    T = 2*pi;                                                   %Initial orbit period
    tspan = 0:dt:T;                                             %Integration time span
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;        %Convergence flag
    maxIter = n;        %Maximum number of iterations   
    iter = 1;           %Initial iteration
    
    %Preallocation 
    ds0 = zeros(2,maxIter);     %Vector containing the initial conditions correction
        
    %Main computation 
    while (GoOn) && (iter < maxIter)
        %Proceed with the integration
        [~, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
        
        %Compute error
        e = [S(end,3); S(end,4)];   %Z and Vx at the crossing must be 0
        
        %Compute the correction 
        F = cr3bp_equations(mu, direction, flagVar, 0, S(end,:).');         %Vector field at T/2
        Phi = reshape(S(end,m+1:end), [m m]);                               %Build the monodromy matrix at T/2
        A = [Phi(3:4,5) Phi(3:4,6)] ...
            -(1/S(end,5)*[S(end,6); F(4)]*[Phi(2,5) Phi(2,6)]);
        ds0(:,iter) = A\e;                                                  %Compute the variation
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            seed(5:6) = seed(5:6)-ds0(:,iter);    %Update initial conditions
            iter = iter+1;                        %Update iteration
        end       
    end
    
    %Integrate the whole trayectory 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Disable crossing event
    tspan = 0:dt:(2*dt*size(S,1));                              %Integrate the orbit for a whole period
    
    [t, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;          %Trajectory
    xf.Period = t(end);         %Orbit period
    
    %Ouput differential correction scheme convergence results
    state = ~GoOn;
end

%Compute periodic orbits using the XZ symmetry
function [xf, state] = SymPlane_scheme(mu, seed, n, tol)
    %Constants 
    m = 6;      %Phase space dimension 
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == 6) || (size(seed,1) == 6)
        %Restrict the seed to the initial conditions
        if (size(seed,2) == 6)
            if (size(seed,1) ~= 1)
                seed = shiftdim(seed(1,:));
            else
                seed = seed.';
            end
        elseif (size(seed,1) == 6) 
            if (size(seed,2) ~= 1)
                seed = shiftdim(seed(:,1));
            end
        end
    else
        disp('No valid initial conditions.');
        xf = []; 
        state = false; 
        return;
    end
    
    %Ensure required initial conditions
    seed(2) = 0;    %Null Y coordinate 
    seed(4) = 0;    %Null Vx 
    seed(6) = 0;    %Null Vz
    
    %Augment initial conditions with the initial STM 
    Phi = eye(m); 
    Phi = reshape(Phi, [m^2 1]); 
    seed(m+1:m+m^2) = Phi;
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, ... 
                     'Events', @(t,s)x_crossing(t,s));          %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    T = 2*pi;                                                   %Initial orbit period
    tspan = 0:dt:T;                                             %Integration time span
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;        %Convergence flag
    maxIter = n;        %Maximum number of iterations   
    iter = 1;           %Initial iteration
    
    %Preallocation 
    ds0 = zeros(2,maxIter);     %Vector containing the initial conditions correction
        
    %Main computation 
    while (GoOn) && (iter < maxIter)
        %Proceed with the integration
        [~, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
        
        %Compute error
        e = [S(end,4); S(end,6)];   %Vx and Vz at the crossing must be 0
        
        %Compute the correction 
        F = cr3bp_equations(mu, direction, flagVar, 0, S(end,:).');         %Vector field at T/2
        Phi = reshape(S(end,m+1:end), [m m]);                               %Build the monodromy matrix at T/2
        A = [Phi(4,3) Phi(4,5); Phi(6,3) Phi(6,5)] ...
            -(1/S(end,5)*[F(4); F(6)]*[Phi(2,3) Phi(2,5)]);
        ds0(:,iter) = A\e;                                                  %Compute the variation
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            seed(3) = seed(3)-ds0(1,iter);        %Update initial conditions
            seed(5) = seed(5)-ds0(2,iter);        %Update initial conditions
            iter = iter+1;                        %Update iteration
        end       
    end
    
    %Integrate the whole trayectory 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Disable crossing event
    tspan = 0:dt:(2*dt*size(S,1));                              %Integrate the orbit for a whole orbit
    [t, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;          %Trajectory
    xf.Period = t(end);         %Orbit period
    
    %Ouput differential correction scheme convergence results
    state = ~GoOn;
end

%Compute periodic orbits using the double X-XZ symmetry
function [xf, state] = SymDouble_scheme(mu, seed, n, tol)
    %Constants 
    m = 6;      %Phase space dimension 
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == 6) || (size(seed,1) == 6)
        %Restrict the seed to the initial conditions
        if (size(seed,2) == 6)
            if (size(seed,1) ~= 1)
                seed = shiftdim(seed(1,:));
            else
                seed = seed.';
            end
        elseif (size(seed,1) == 6) 
            if (size(seed,2) ~= 1)
                seed = shiftdim(seed(:,1));
            end
        end
    else
        disp('No valid initial conditions.');
        xf = []; 
        state = false; 
        return;
    end
    
    %Ensure required initial conditions
    seed(2) = 0;    %Null Y coordinate 
    seed(3) = 0;    %Null Z coordinate
    seed(4) = 0;    %Null Vx 
    
    %Augment initial conditions with the initial STM 
    Phi = eye(m); 
    Phi = reshape(Phi, [m^2 1]); 
    seed(m+1:m+m^2) = Phi;
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, ... 
                     'Events', @(t,s)x_crossing(t,s));          %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    T = 2*pi;                                                   %Initial orbit period
    tspan = 0:dt:T;                                             %Integration time span
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;        %Convergence flag
    maxIter = n;        %Maximum number of iterations   
    iter = 1;           %Initial iteration
    
    %Preallocation 
    ds0 = zeros(2,maxIter);     %Vector containing the initial conditions correction
        
    %Main computation 
    while (GoOn) && (iter < maxIter)
        %Proceed with the integration
        [~, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
        
        %Compute error
        e = [S(end,4); S(end,6)];   %Vx and Vz at the crossing must be 0
        
        %Compute the correction 
        F = cr3bp_equations(mu, direction, flagVar, 0, S(end,:).');   %Vector field at T/4
        Phi = reshape(S(end,m+1:end), [m m]);                         %Build the monodromy matrix at T/4
        A(1,1) = Phi(4,5)-(F(4)/S(end,5))*Phi(2,5);
        A(1,2) = Phi(4,6)-(F(4)/S(end,5))*Phi(2,6);
        A(2,1) = Phi(6,5)-(F(6)/S(end,5))*Phi(2,5);
        A(2,2) = Phi(6,6)-(F(6)/S(end,5))*Phi(2,6);
        ds0(:,iter) = A\e;                                            %Compute the variation
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            seed(5:6) = seed(5:6)-ds0(:,iter);    %Update initial conditions
            iter = iter+1;                        %Update iteration
        end       
    end
    
    %Integrate the whole trayectory 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Disable crossing event
    tspan = 0:dt:(4*dt*size(S,1));                              %Integrate the orbit for a whole period
    [t, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;          %Trajectory
    xf.Period = t(end);         %Orbit period
    
    %Ouput differential correction scheme convergence results
    state = ~GoOn;
end

%Compute periodic orbits using the symmetry XZ-XY planes
function [xf, state] = SymDoublePlane_scheme(mu, seed, n, tol)
    %Constants 
    m = 6;      %Phase space dimension 
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == 6) || (size(seed,1) == 6)
        %Restrict the seed to the initial conditions
        if (size(seed,2) == 6)
            if (size(seed,1) ~= 1)
                seed = shiftdim(seed(1,:));
            else
                seed = seed.';
            end
        elseif (size(seed,1) == 6) 
            if (size(seed,2) ~= 1)
                seed = shiftdim(seed(:,1));
            end
        end
    else
        disp('No valid initial conditions.');
        xf = []; 
        state = false; 
        return;
    end
    
    %Ensure required initial conditions
    seed(2) = 0;    %Null Y coordinate 
    seed(4) = 0;    %Null Vx 
    seed(6) = 0;    %Null Vz
    
    %Augment initial conditions with the initial STM 
    Phi = eye(m); 
    Phi = reshape(Phi, [m^2 1]); 
    seed(m+1:m+m^2) = Phi;
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, ... 
                     'Events', @(t,s)x_crossing(t,s));          %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    T = 2*pi;                                                   %Initial orbit period
    tspan = 0:dt:T;                                             %Integration time span
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;        %Convergence flag
    maxIter = n;        %Maximum number of iterations   
    iter = 1;           %Initial iteration
    
    %Preallocation 
    ds0 = zeros(2,maxIter);     %Vector containing the initial conditions correction
        
    %Main computation 
    while (GoOn) && (iter < maxIter)
        %Proceed with the integration
        [~, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
        
        %Compute error
        e = [S(end,4); S(end,6)];   %Vx and Vz at the crossing must be 0
        
        %Compute the correction 
        F = cr3bp_equations(mu, direction, flagVar, 0, S(end,:).');         %Vector field at T/2
        Phi = reshape(S(end,m+1:end), [m m]);                               %Build the monodromy matrix at T/2
        A = [Phi(4,1) Phi(4,6); Phi(6,1) Phi(6,6)] ...
            -(1/S(end,5)*[F(4); F(6)]*[Phi(2,1) Phi(2,6)]);
        ds0(:,iter) = A\e;                                                  %Compute the variation
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            seed(1) = seed(1)-ds0(1,iter);        %Update initial conditions
            seed(6) = seed(6)-ds0(2,iter);        %Update initial conditions
            iter = iter+1;                        %Update iteration
        end       
    end
    
    %Integrate the whole trayectory 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Disable crossing event
    tspan = 0:dt:(2*dt*size(S,1));                              %Integrate the orbit for a whole orbit
    [t, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;          %Trajectory
    xf.Period = t(end);         %Orbit period
    
    %Ouput differential correction scheme convergence results
    state = ~GoOn;
end

%Compute planar periodic orbits -for Lyapunov orbits-
function [xf, state] = SymPlanar_scheme(mu, seed, n, tol)
    %Constants 
    m = 6;      %Phase space dimension 
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == 6) || (size(seed,1) == 6)
        %Restrict the seed to the initial conditions
        if (size(seed,2) == 6)
            if (size(seed,1) ~= 1)
                seed = shiftdim(seed(1,:));
            else
                seed = seed.';
            end
        elseif (size(seed,1) == 6) 
            if (size(seed,2) ~= 1)
                seed = shiftdim(seed(:,1));
            end
        end
    else
        disp('No valid initial conditions.');
        xf = []; 
        state = false; 
        return;
    end
    
    %Ensure motion on the XY synodic plane
    seed(2) = 0;    %Null Y coordinate 
    seed(3) = 0;    %Null Z coordinate
    seed(4) = 0;    %Null Vx 
    seed(6) = 0;    %Null Vz 
    
    %Augment initial conditions with the initial STM 
    Phi = eye(m); 
    Phi = reshape(Phi, [m^2 1]); 
    seed(m+1:m+m^2) = Phi;
    
    %Set up integration 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22, ... 
                     'Events', @(t,s)x_crossing(t,s));          %Integration conditions and tolerances                     
    dt = 1e-3;                                                  %Time step
    T = 2*pi;                                                   %Initial orbit period
    tspan = 0:dt:T;                                             %Integration time span
    direction = 1;                                              %Forward integration
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;        %Convergence flag
    maxIter = n;        %Maximum number of iterations   
    iter = 1;           %Initial iteration
    
    %Preallocation 
    ds0 = zeros(1,maxIter);     %Vector containing the initial conditions correction
        
    %Main computation 
    while (GoOn) && (iter < maxIter)
        %Proceed with the integration
        [~, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
        
        %Compute error
        e = S(end,4);   %Vx at the crossing must be 0
        
        %Compute correction 
        F = cr3bp_equations(mu, direction, flagVar, 0, S(end,:).');   %Vector field at T/2
        Phi = reshape(S(end,m+1:end), [m m]);                         %Build the monodromy matrix at T/2
        ds0(iter) = e/(Phi(4,5)-(F(4)/S(end,5))*Phi(2,5));            %Compute the variation
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            seed(5) = seed(5)-ds0(iter);    %Update initial conditions
            iter = iter+1;                  %Update iteration
        end       
    end
    
    %Integrate the whole trayectory 
    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);      %Disable crossing event
    tspan = 0:dt:(2*dt*size(S,1));
    [t, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;          %Trajectory
    xf.Period = t(end);         %Orbit period
    
    %Ouput differential correction scheme convergence results
    state = ~GoOn;
end

%Compute periodic orbits using multiple shooting and energy-continuity constraint
function [xf, state] = MSPeriodic_scheme(mu, seed, n, tol, varargin)
    %Assign undeclared local inputs if any. Sanity check 
    if (isempty(varargin{1}))
       disp('No valid inputs. Correction is about to finish.');
       xf = []; 
       state = false;
       return;
    else
        local_inputs = varargin{1};
        nodes = local_inputs{1};            %Nodes to compute
        T = local_inputs{2};                %Initial period of the orbit
        
        if (nodes < 2) 
            disp('No valid inputs. Correction is about to finish.'); 
            xf = []; 
            state = false;
            return;
        end
    end
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == 6) || (size(seed,1) == 6)
        if (size(seed,2) == 6)
            seed = seed.';          %Accomodate new format
        end
    else
        disp('No valid initial conditions.');
        xf = []; 
        state = false; 
        return;
    end
    
    %Constants 
    m = 6;                          %Phase space dimension 
    Phi = eye(m);                   %Initial STM  
    Phi = reshape(Phi, [m^2 1]);    %Initial STM 
    dt = 1e-4;                      %Integration time step
    h = fix(size(seed,2)/nodes)-1;  %Temporal index step
    Dt = T/nodes;                   %Time step
    constraints = 6;                %Additional constraints to continuity
        
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
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;                                                %Convergence flag
    maxIter = n;                                                %Maximum number of iterations   
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
            S0 = [shiftdim(internalSeed(m*(i-1)+1:m*i)); Phi];
            [~, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, S0, options);
            F = cr3bp_equations(mu, direction, flagVar, 0, S(end,:).');          %Vector field
            
            %Build the covariance matrix                                         %Vector field matrix
            if (i ~= nodes)
                %Continuity constraint
                A(m*(i-1)+1:m*i,m*(i-1)+1:m*i) = reshape(S(end,m+1:end),[m m]);  %Subarc STM
                A(m*(i-1)+1:m*i,m*i+1:m*(i+1)) = -eye(m);                        %Continuity constraint matrix
                B(m*(i-1)+1:m*i,i) = F(1:m);                                     %Dynamics matrix
            else
                %Periodicity constraint
                STM = reshape(S(end,m+1:end),[m, m]);                            %Subarc STM
                A(end-m+1:end-1,end-m+1:end) = [STM(1:4,:); STM(1,:)];           %Constraint matrix
                A(end-m+1:end-1,1:m) = -[eye(4) zeros(4,2); zeros(1,5) 1];       %Constraint matrix          
                %Jacobi Constant constraint
                A(end,end-m+1:end) = -jacobi_gradient(mu, S(end,1:m).').';       %Constraint matrix
                A(end,1:m) = -jacobi_gradient(mu, internalSeed(1:m)).';          %Constraint matrix
            end     
            
            %Compute the error
            if (i ~= nodes)
                e(m*(i-1)+1:m*i) = shiftdim(S(end,1:m).'-internalSeed(m*i+1:m*(i+1)));  %Continuity constraint
            else
                dR = shiftdim(S(end,1:m).'-internalSeed(1:m));
                e(end-m+1:end-1) = [dR(1:4); dR(6)];                                                %Periodicity constraint
                e(end) = jacobi_constant(mu, internalSeed(1:m))-jacobi_constant(mu, S(end,1:m).');  %Jacobi Constant constraint
            end
        end
        
        %Full covariance matrix 
        C = [A B];
                
        %Compute the correction 
        ds0(:,iter) = C.'*(C*C.')^(-1)*e;               %Compute the variation (under-determined case)
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            internalSeed = internalSeed-ds0(:,iter);    %Update initial conditions
            iter = iter+1;                              %Update iteration
        end       
    end
    
    %Integrate the whole trayectory
    tspan = 0:dt:sum(internalSeed(end-nodes+1:end))+Dt;
    seed = [shiftdim(internalSeed(1:m)); Phi];                  
    [t, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;                           %Trajectory
    xf.Period = t(end);                          %Orbit period
        
    %Ouput differential correction scheme convergence results
    state = ~GoOn;
end

%Compute periodic orbits using multiple shooting and fixed Jacobi Constant value 
function [xf, state] = MSJacobi_scheme(mu, seed, n, tol, varargin)
    %Assign undeclared local inputs if any. Sanity check 
    if (isempty(varargin{1}))
       disp('No valid inputs. Correction is about to finish.');
       xf = []; 
       state = false;
       return;
    else
        local_inputs = varargin{1};
        nodes = local_inputs{1};            %Nodes to compute
        T = local_inputs{2};                %Initial period of the orbit
        Cref = local_inputs{3};             %Jacobi constant reference value
        
        if (nodes < 2) 
            disp('No valid inputs. Correction is about to finish.'); 
            xf = []; 
            state = false;
            return;
        end
    end
    
    %Sanity check on initial conditions dimension
    if (size(seed,2) == 6) || (size(seed,1) == 6)
        if (size(seed,2) == 6)
            seed = seed.';          %Accomodate new format
        end
    else
        disp('No valid initial conditions.');
        xf = []; 
        state = false; 
        return;
    end
    
    %Constants 
    m = 6;                          %Phase space dimension 
    Phi = eye(m);                   %Initial STM  
    Phi = reshape(Phi, [m^2 1]);    %Initial STM 
    dt = 1e-4;                      %Integration time step
    h = fix(size(seed,2)/nodes)-1;  %Temporal index step
    Dt = T/nodes;                   %Time step
    constraints = 6;                %Additional constraints to continuity
        
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
    flagVar = true;                                             %Integrate variational equations
    
    %Set up differential correction scheme
    GoOn = true;                                                %Convergence flag
    maxIter = n;                                                %Maximum number of iterations   
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
            S0 = [shiftdim(internalSeed(m*(i-1)+1:m*i)); Phi];
            [~, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, S0, options);
            F = cr3bp_equations(mu, direction, flagVar, 0, S(end,:).');          %Vector field
            
            %Build the covariance matrix                                         %Vector field matrix
            if (i ~= nodes)
                %Continuity constraint
                A(m*(i-1)+1:m*i,m*(i-1)+1:m*i) = reshape(S(end,m+1:end),[m m]);  %Subarc STM
                A(m*(i-1)+1:m*i,m*i+1:m*(i+1)) = -eye(m);                        %Continuity constraint matrix
                B(m*(i-1)+1:m*i,i) = F(1:m);                                     %Dynamics matrix
            else
                %Periodicity constraint
                STM = reshape(S(end,m+1:end),[m, m]);                            %Subarc STM
                A(end-m+1:end-1,end-m+1:end) = [STM(1:4,:); STM(1,:)];           %Constraint matrix
                A(end-m+1:end-1,1:m) = -[eye(4) zeros(4,2); zeros(1,5) 1];       %Constraint matrix          
                %Jacobi Constant constraint
                A(end,end-m+1:end) = -jacobi_gradient(mu, S(end,1:m).').';       %Constraint matrix
            end     
            
            %Compute the error
            if (i ~= nodes)
                e(m*(i-1)+1:m*i) = shiftdim(S(end,1:m).'-internalSeed(m*i+1:m*(i+1)));  %Continuity constraint
            else
                dR = shiftdim(S(end,1:m).'-internalSeed(1:m));
                e(end-m+1:end-1) = [dR(1:4); dR(6)];                                    %Periodicity constraint
                e(end) = Cref-jacobi_constant(mu, S(end,1:m).');                        %Jacobi Constant constraint
            end
        end
        
        %Full covariance matrix 
        C = [A B];
                
        %Compute the correction 
        ds0(:,iter) = C.'*(C*C.')^(-1)*e;               %Compute the variation (under-determined case)
        
        %Convergence analysis 
        if (norm(e) <= tol)
            GoOn = false;
        else
            internalSeed = internalSeed-ds0(:,iter);    %Update initial conditions
            iter = iter+1;                              %Update iteration
        end       
    end
    
    %Integrate the whole trayectory
    tspan = 0:dt:sum(internalSeed(end-nodes+1:end))+Dt;
    seed = [shiftdim(internalSeed(1:m)); Phi];                  
    [t, S] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), tspan, seed, options);
    
    %Ouput corrected trajectory 
    xf.Trajectory = S;                           %Trajectory
    xf.Period = t(end);                          %Orbit period
        
    %Ouput differential correction scheme convergence results
    state = ~GoOn;
end
