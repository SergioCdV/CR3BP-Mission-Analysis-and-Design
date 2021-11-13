%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/1/21
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

function [u] = MLQR_control(mu, Sg, Sn, Q, M)
    %Approximation 
    n = 6;                                          %Dimension of the state vector
    
    %Preallocation 
    u = zeros(3,size(Sn,1));                        %Control law      

    for i = 1:size(Sn,1)
        %Compute the actual state 
        X(1) = jacobi_constant(mu, Sn(1,1:n).');    %Energy level of the spacecraft
        X(2) = unstable_component(Sn(1,:).');       %Unstable coordinate in the Floquet basis

        %Numerically compute the derivative of the unstable Floquet mode
        dalpha = numerical_jacobian(1, @(Sn)unstable_component(Sn), Sn(i,:).');
        dalpha = dalpha(1:n);

        %Compute the state matrix
        Monodromy = reshape(Sn(1,n+1:end), [n n]);
        [E, lambda] = eig(Monodromy);               %Eigenspectrum of the STM 
        A = [0 0; 0 lambda(1,1)];                   %State transition matrix
        J = abs_jacobian(mu, Sn(1,1:n).');
        J = E*J*E^(-1); 

        %Compute the control matrix 
        B = [zeros(3,3); eye(3)];

        %Compute the gradient of the integrals 
        C = jacobi_gradient(mu, Sn(1,1:n).');
        V = [C.'; dalpha]*B;

        %Compute the feedback control law
        if (abs(lambda-1) <= 1)
            u(:,i) = zeros(3,1); 
        else
            [K,~,~] = lqr(A,V,Q,M);
            u(:,i) = -K*(X-Sg(i,:)).';
        end
    end
end

%% Auxiliary functions 
%Compute the state space derivative of the unstable component of the state vector
function [dalpha] = unstable_component(Sn)
    %Constants 
    n = 6;              %State dimension
        
    %Compute the gradient of the STM and its eigenvalues with respect to the state variable
    Monodromy = reshape(Sn(n+1:end), [n n]);
    [E, ~] = eig(Monodromy);                            %Eigenspectrum of the STM 
    for j = 1:size(E,2)
        E(:,j) = E(:,j)/norm(E(:,j));                   %Compute the Floquet Modes
    end
    chi = E*Sn(1:n);                                    %Linearized coordinates in the Floquet basis
    dalpha = dot(E(:,1), chi);                          %Projection onto the unstable manifold
end