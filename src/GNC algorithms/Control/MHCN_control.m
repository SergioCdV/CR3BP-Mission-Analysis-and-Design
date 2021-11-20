%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 20/11/21
% File: MHCN_control.m 
% Issue: 0 
% Validated: 

%% Manifold Hybrid Control %%
% This script contains the function to compute the control law by means of a hybrid, manifold-based provide low-thrust stationkeeping.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - array Sg, the guidance law to follow (reference orbit Jacobi
%           Constant)
%         - array Sn, the system state (with the STM)
%         - matrices Q and M, penalizing on the state error and the control
%           effort

% Output: - vector u, the computed control law

% New versions: 

function [u] = MHCN_control(mu, t, T, Sg, Sn, Q, M)
    %Approximation 
    n = 6;                                          %Dimension of the state vector
    
    %Preallocation 
    u = zeros(3,size(Sn,1));                        %Control law   

    %Compute the error state 
    dSn = zeros(size(Sn));
    dSn(:,1:n) = Sn(:,1:n)-Sg(:,1:n);
    dSn(:,n+1:end) = Sn(:,n+1:end);

    %Compute the control matrix 
    B = [zeros(3,3); eye(3)];

    for i = 1:size(Sn,1)
        %Compute the actual state 
        X(1) = unstable_component(t(i), T, dSn(i,:).');    %Unstable coordinate in the Floquet basis
        X(2) = stable_component(t(i), T, dSn(i,:).');      %Stable coordinate in the Floquet basis

        %Compute the feedback control law
        if (X(1) > 1e-3)
            %Eigenspectrum of the STM
            Monodromy = reshape(Sn(i,n+1:end), [n n]);
            [E, lambda] = eig(Monodromy);  

            %Numerically compute the derivative of the unstable Floquet mode
            dalpha = exp(-t(i)/T*log(lambda(1,1)))*(E(:,1).'+dSn(i,1:n));      %Derivative of the projection onto the unstable manifold;
    
            %Compute the state matrix
            A = [0 0; 0 lambda(1,1)];                   %State transition matrix
    
            %Compute the gradient of the integrals 
            C = jacobi_gradient(mu, Sn(i,1:n).');
            V = [C.'; dalpha]*B;
            [K,~,~] = lqr(A,V,Q,M);
            u(:,i) = -K*(X-Sg(i,end-1:end)).';
        else
            u(:,i) = zeros(3,1);
        end
    end
end

%% Auxiliary functions 
%Compute the unstable component of the state vector
function [alpha] = unstable_component(t, T, Sn)
    %Constants 
    n = 6;              %State dimension
        
    %Compute the gradient of the STM and its eigenvalues with respect to the state variable
    Monodromy = reshape(Sn(n+1:end), [n n]);
    [E, lambda] = eig(Monodromy);                                     %Eigenspectrum of the STM 
    alpha = exp(-t/T*log(lambda(1,1)))*dot(E(:,1),Sn(1:n));           %Projection onto the unstable manifold
end

%Compute the stable component of the state vector
function [alpha] = stable_component(t, T, Sn)
    %Constants 
    n = 6;              %State dimension
        
    %Compute the gradient of the STM and its eigenvalues with respect to the state variable
    Monodromy = reshape(Sn(n+1:end), [n n]);
    [E, lambda] = eig(Monodromy);                                     %Eigenspectrum of the STM 
    alpha = exp(-t/T*log(lambda(2,2)))*dot(E(:,2),Sn(1:n));           %Projection onto the stable manifold
end

%Compute the state space derivative of the unstable component of the state vector
function [dalpha] = dunstable_component(mu, t, T, Sn)
    %Constants 
    n = 6;              %State dimension
        
    %Compute the gradient of the STM and its eigenvalues with respect to the state variable
    Monodromy = reshape(Sn(n+1:end), [n n]);
    [E, lambda] = eig(Monodromy);                                %Eigenspectrum of the STM 
    J = abs_jacobian(mu, Sn);
    dalpha = exp(-t/T*log(lambda(1,1)))*(E(:,1)).';      %Derivative of the projection onto the unstable manifold
end