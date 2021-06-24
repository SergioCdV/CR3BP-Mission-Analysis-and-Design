%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/06/21
% File: libration_curves.m 
% Issue: 0 
% Validated: 

%% Libration curves %%
% For a given gravitational parameter mu and target motion, this function computes the 
% relative equilibrium curves.

% Inputs: - scalar mu, the reduced gravitational parameter of the system.
%         - array rt, the time evolution of the target motion

% Outputs: - array L, containing in one matrix the position of the five libration 
%            points of the system (with each column being a position vector) in L1, L2, 
%            L3, L4 and L5 order, as well as the collinear gamma distance for the first three.

% New versions: 

function [L] = libration_curves(mu, rt)
    %Preallocation 
    L = zeros(size(rt,1), 15);               %Equilibrium curves preallocation
    
    %Initial guess 
    absL = libration_points(mu);             %Absolute libration points
    
    %Compute the curves for each initial guess
    for i = 1:size(rt,1)
        for j = 1:size(absL,2)
            L(i,1+3*(j-1):3*j) = curve_solver(mu, rt(i,1:3).', absL(1:3,j));
        end
    end
end

%% Auxiliary functions 
%Solve the equilibriuum gradient
function [Lc] = curve_solver(mu, rt, L)
    %Constants 
    R(1:3,1) = [-mu; 0; 0];        %Location of the absolute first primary
    R(1:3,2) = [1-mu; 0; 0];       %Location of the absolute second primary
    Rr(1:3,1) = R(1:3,1)-rt;       %Location of the unsteady first primary
    Rr(1:3,2) = R(1:3,2)-rt;       %Location of the unsteady second primary
    
    mu_r(1) = 1-mu;                %Gravitational parameter of the first primary
    mu_r(2) = mu;                  %Gravitational parameter of the second primary
    
    %Newton solver setup
    tol = 1e-3;                    %Convergence tolerance
    GoOn = true;                   %Convergence flag
    maxIter = 100;                 %Maximum number of iterations
    iter = 1;                      %Initial iterations
    
    %Preallocation 
    gradient = zeros(3,1);         %Equilibrium condition
    Lc = zeros(3,maxIter);         %Preallocation of the solution
    Lc(:,iter) = [L(1:2); rt(3)];  %Initial guess        
    
    %Newton's solver
    while (GoOn) && (iter < maxIter)
        %Configuration space variables 
        x = Lc(1,iter);            %Relative synodic x coordinate
        y = Lc(2,iter);            %Relative synodic y coordinate
        z = Lc(3,iter);            %Relative synodic z coordinate
        
        %Function and Jacobian evaluation
        f = zeros(3,1);
        for i = 1:length(mu_r)
            gradient(1) = mu_r(i)*(-Rr(1,i)/norm(Rr(:,i))^3-(x-Rr(1,i))/norm(Lc(:,iter)-Rr(:,i))^3);
            gradient(2) = mu_r(i)*(-Rr(2,i)/norm(Rr(:,i))^3-(y-Rr(2,i))/norm(Lc(:,iter)-Rr(:,i))^3);
            gradient(3) = mu_r(i)*(-Rr(3,i)/norm(Rr(:,i))^3-(z-Rr(3,i))/norm(Lc(:,iter)-Rr(:,i))^3);
            f = f + gradient;
        end
        f = f + [x; y; 0];
        
        Phi = rel_jacobian(mu, [rt; zeros(3,1); Lc(:,iter); zeros(3,1)]); 
        Phi = Phi(4:6,1:3);
        
        %Newton residual 
        ds = -Phi\f;
            
        %Update the equilibrium point
        Lc(:,iter+1) = Lc(:,iter) + ds; 
        
        %Convergence flag
        if (norm(ds) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end
    
    %Select the last converged value 
    Lc = Lc(:,iter);
end