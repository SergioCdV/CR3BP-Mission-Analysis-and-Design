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

    %Model coefficients 
    mup(1) = 1-mu;                                  %Reduced gravitational parameter of the first primary 
    mup(2) = mu;                                    %Reduced gravitational parameter of the second primary 
    R(:,1) = [-mu; 0; 0];                           %Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];                          %Synodic position of the second primary

    %Linear model matrices
    B = [zeros(n/2); eye(n/2)];                     %Linear model input matrix 
    
    %Preallocation 
    u = zeros(3,size(Sn,1));                        %Control law      

    %Linear state model
    A = zeros(2,6);  
    
    for i = 1:size(Sn,1)
        %Compute the gradient of the STM and its eigenvalues with respect to the state variable

        %Compute the gradient of the integrals 
        C = jacobi_gradient(mu, Sn(1,1:n).');
        V = [C; dalpha]*B;

        %Compute the feedback control law
        [K,~,~] = lqr(A,V,Q,M);
        u(:,i) = -K*(Sn(i,:)-Sg(i,:)).';
    end
end