%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: invariant_manifold.m 
% Issue: 0 
% Validated: 

%% Invariant manifold %%
% This function contains the algorithm to compute the invariant manifolds associated with 
% a particular solution orbit.

% Inputs: - scalar mu, the reduced gravitational parameter of the system.
%         - string manifold, selecting the manifold to compute, either 'U' for the unstable one or 
%         - 'S', for the stable one. 
%         - character branch, 'R' for the right branch or 'L' for the left
%           one.
%         - vector r, a temporal evolution of the orbit in an MxN matrix,
%           with N = 42, containing the phase space vector and the STM
%           evolutions, and M being the number of temporal steps to
%           complete a orbital period.
%         - scalar rho, a number of fibers/trajectories to compute on the
%           manifold.
%         - scalar time, to integrate the dynamics

% Output: - vector field M, containing the manifold evolution for each
%           fiber. 

% Methods: integration of the dynamics and use of the monodromy matrix and
%          its propagation all along the orbit.

% New versions:

function [M] = invariant_manifold(mu, manifold, branch, r, rho, tspan)
    %General constants 
    flagVar = false;             %No STM integration needed
    n = 6;                       %Phase space dimension
    epsilon = 1e-5;              %Displacement of the initial conditions  
    T = size(r,1);               %Orbit period in nondimensinal units
    
    %Integration tolerances
    RelTol = 2.25e-14; 
    AbsTol = 1e-22;
    options = odeset('RelTol', RelTol, 'AbsTol', AbsTol);
    
    %Integration and manifold parameters     
    if (manifold == 'U')
        direction = 1;          %Forward integration
        eigenV = 1;             %Unstable eigenvector
    elseif (manifold == 'S')
        direction = -1;         %Backward integration
        eigenV = 2;             %Stable eigenvector
    else
        disp('No valid manifold was selected. Try again.'); 
        disp(' ');
        M = [];
    end
    
    if (branch == 'L')
        epsilon = -epsilon;
    end
    
    %Monodromy matrix from the trayectory r
    monodromy = reshape(r(end,n+1:end), n, n);
    
    %Eigenvalues and eigenvectors of the monodromy matrix
    [W, ~] = eig(monodromy);                    %Eigenvector and eigenvalues
    mV0 = W(:,eigenV);                          %Selected eigenvector
    
    %Select insertion points along the orbit from which compute the manifold fibers
    M0 = zeros(rho, n);                         %Preallocate for speed
    h = round((T-1)/rho);                       %Spatial step
    for i = 1:rho
        %Position and time in the periodic orbit 
        orbitT = (i-1)*h+1;                     %Orbit independent variables (t or theta, using a time law)
        orbitX0 = r(orbitT,1:n);                %Orbit point position
        
        %STM at that point and time
        Phi = reshape(r(orbitT,n+1:end), [n n]);
        
        %Initial conditions of the manifold fiber
        mV = Phi*mV0;                           %Pushed and propagate the selected eigenvector
        mV = mV/norm(mV);                       %Normalized propagated selected eigenvector
        M0(i,:) = orbitX0+epsilon*mV.';         %Manifold initial condition
    end
    
    %Complete manifold integration 
    M = zeros(rho, length(tspan), size(M0,2));
    for i = 1:rho
    	[~, auxM] = ode113(@(t,s)cr3bp_equations(mu, direction, flagVar, t, s), ...
                            tspan, M0(i,:), options);
        M(i,:,:) = auxM;
    end 
end