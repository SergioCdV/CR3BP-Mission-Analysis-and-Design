%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/11/21
% File: MLQR_control.m 
% Issue: 0 
% Validated: 

%% Manifold Linear Quadratic Regulator Control %%
% This script contains the function to compute the control law by means of an LQR controller to provide low-thrust stationkeeping.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - array Sg, the guidance law to follow (reference orbit Jacobi
%           Constant)
%         - array Sn, the system state (with the STM)
%         - matrices Q and M, penalizing on the state error and the control
%           effort

% Output: - vector u, the computed control law

% New versions: 

function [u] = MLQR_control(mu, t, T, Sg, Sn, Q, M)
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
        X(1) = jacobi_constant(mu, Sn(i,1:n).');           %Energy level of the spacecraft
        X(2) = unstable_component(t(i), T, dSn(i,:).');    %Unstable coordinate in the Floquet basis
        X(3) = stable_component(t(i), T, dSn(i,:).');      %Stable coordinate in the Floquet basis

        %Compute the feedback control law
        if (mod(t(i),T) > 0.9 || mod(t(i),T) == 0)
            u(:,i) = zeros(3,1); 
        else
            %Eigenspectrum of the STM
            Monodromy = reshape(Sn(i,n+1:end), [n n]);
            J = abs_jacobian(mu,Sn(i,1:n));
            [E, lambda] = eig(Monodromy);  

            %Numerically compute the derivative of the unstable Floquet mode
            dtE = J*E-E*lambda;
            daE = dtE(:,1); 
            dE = daE./Sn(i,1:n);
            dalpha = exp(-t(i)/T*log(lambda(1,1)))*(E(:,1).'+dSn(i,1:n)*dE);      %Derivative of the projection onto the unstable manifold
            dbE = dtE(:,2); 
            dE = dbE./Sn(i,1:n);
            dbeta = exp(-t(i)/T*log(lambda(1,1)))*(E(:,2).'+dSn(i,1:n)*dE);       %Derivative of the projection onto the stable manifold
    
            %Compute the state matrix
            A = [0 0 0; 0 log(lambda(1,1))/T 0; 0 0 log(lambda(2,2))/T];          %State transition matrix
    
            %Compute the gradient of the integrals 
            C = jacobi_gradient(mu, Sn(i,1:n).');
            V = [C.'; dalpha; dbeta]*B;
            [K,~,~] = lqr(A,V,Q,M);
            e = (X-Sg(i,end-2:end)).';
            u(:,i) = -K*e;
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