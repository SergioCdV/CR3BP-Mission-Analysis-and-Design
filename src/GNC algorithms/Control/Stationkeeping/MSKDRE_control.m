%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/11/21
% File: MSKDRE_control.m 
% Issue: 0 
% Validated: 

%% Manifold Stationkeeping State Dependent Ricatti Equation %%
% This script contains the function to compute the control law by means of an SDRE controller to provide low-thrust stationkeeping.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar t, the current simulation time
%         - matrix F, the initial conditions of the periodic orbit Floquet
%           basis
%         - array St, the system target state evolution
%         - scalar Jref, the referecen Jacobi Constant value
%         - array Sn, the system state (with the STM)
%         - matrices Q and R, to compute the LQR controller

% New versions: 

function [u] = MSKDRE_control(mu, t, T, L, F, St, Sg, Sn, Q, M)
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
        J = jacobi_constant(mu, Sr.');          %Absolute Jacobi Constant 
        e(2) = J-Sg;                            %Jacobi Constant error
        e = e.';

        %Compute the gradient of the integrals 
        V = invE*B;                             %Control input matrix
        V = V(1,:);                             %Control input matrix
        dJ = jacobi_gradient(mu,Sr.').';        %Gradient of the Jacobi Constant 
        V(2,:) = dJ*B;                          %Gradient of the Jacobi Constant with respect to the Floquet variables
        C = dJ*E*L;                             %Gradient of the Jacobi Constant with respect to the relative Floquet variables 
        A(2,1) = C(1);                          %Gradient of the Jacobi Constant with respect to the unstable relative Floquet variable

        %SDRE design
        C = ctrb(A,V);                          %Controlabillity matrix

        %Control law
        if (rank(C) ~= max(size(A)))
            u(:,i) = zeros(3,1);           
        else
            [K,~,~] = lqr(A,V,Q,M);  
            u(:,i) = -real(K*e);           
        end
    end
end