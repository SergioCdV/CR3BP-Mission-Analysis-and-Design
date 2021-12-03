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
%         - array Sn, the system state (with the STM)
%         - matrices Q and M, penalizing on the state error and the control
%           effort

% Output: - vector u, the computed control law

% New versions: 

function [u] = MLQR_control(mu, t, T, L, F, St, Sg, Sn, Q, M)
    %Approximation 
    n = 6;                                          %Dimension of the state vector
    
    %Preallocation 
    u = zeros(3,size(Sn,1));                        %Control law  

    %Compute the control matrix 
    B = [zeros(3,3); eye(3)];

    %Compute the state matrix
    A = [L(1,1) 0; 0 0];

    for i = 1:size(Sn,1)
        %Compute the actual state
        Monodromy = reshape(Sn(i,n+1:n+n^2), [n n]);

        E = Monodromy*F*expm(-L*t(i));          %Floquet basis                                                                                                           
        invE = E^(-1);                          %Floquet change of basis
        e = invE*Sn(i,1:n).';                   %State error
        e = e(1);                               %Unstable component of the error

        %Jacobi Constant constraint
        Sr = St(i,1:n)+Sn(i,1:n);               %Absolute trajectory
        e(2) = jacobi_constant(mu, Sr.');       %Absolute Jacobi Constant 
        e(2) = e(2)-Sg;                         %Jacobi Constant error
        e = e.'; 

        %Compute the gradient of the integrals 
        V = invE*B;                             %Control input matrix
        V = V(1,:);                             %Control input matrix
        dJ = jacobi_gradient(mu,Sr.').';        %Gradient of the Jacobi Constant 
        V(2,:) = dJ*invE*B;                     %Gradient of the Jacobi Constant with respect to the Floquet variables

        %SDRE design
        C = ctrb(A,V);                          %Controlabillity matrix
        [K,~,~] = lqr(A,V,Q,M);  

       % norm(e)

        %Control law
        if (rank(C) ~= max(size(A)))
            u(:,i) = zeros(3,1);           
        else
            u(:,i) = real(-K*e);           
        end
    end
end