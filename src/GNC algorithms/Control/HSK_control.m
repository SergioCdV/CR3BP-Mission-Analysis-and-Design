%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/08/21
% File: HSK_control.m 
% Issue: 0 
% Validated: 08/08/21

%% Hamiltonian Station-keeping Control %%
% This script contains the function to compute the control law by means of an HSK controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - array Sg, the guidance law to follow
%         - array Sn, the system state
%         - matrices Q and M, penalizing on the state error and the control
%           effort

% Output: - vector u, the computed control law

% New versions: 

function [u] = HSK_control(mu, Sg, Sn, Q, R)
    %Approximation 
    n = 6;                                          %Dimension of the state vector

    %Model coefficients 
    B = [zeros(n/2,1); ones(n/2,1)];                %Linear model input matrix 
    
    %Preallocation 
    u = zeros(size(B,1), size(Sg,1));
    
    for i = 1:size(Sg,1)
        %Compute the Jacobi integral associated to the guidance law 
        Jref = jacobi_constant(mu, Sg(i,:).')+0.1;

        %Compute the control law 
        J = jacobi_constant(mu, Sn(i,:).');             %Jacobi Constant of the instantenous state
        dJ = jacobi_gradient(mu, Sn(i,:).');            %Gradient of the Jacobi Constant of the instantenous state 
        Bj = dot(dJ,B);                                 %Projection of the gradient of the Jacobi Constant on the state space
        U = -sign(Bj)*sqrt(Q/R)*(J-Jref);               %Scalar control law 
        u(:,i) = B*U;                                   %Full control vector     
    end
end