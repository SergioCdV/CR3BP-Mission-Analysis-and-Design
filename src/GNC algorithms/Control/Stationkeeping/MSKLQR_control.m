%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/11/21
% File: MLQR_control.m 
% Issue: 0 
% Validated: 

%% Manifold Stationkeeping Linear Quadratic Regulator Control %%
% This script contains the function to compute the control law by means of an LQR controller to provide low-thrust stationkeeping.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - scalar t, the current simulation time
%         - matrix F, the initial conditions of the periodic orbit Floquet
%           basis
%         - array St, the system target state evolution
%         - scalar Jref, the referecen Jacobi Constant value
%         - array Sn, the system state (with the STM)
%         - matrix K, the system controller

% Output: - vector u, the computed control law

% New versions: 

function [u] = MSKLQR_control(mu, t, T, L, F, St, Jref, Sn, K)
    %Approximation 
    n = 6;                                          %Dimension of the state vector
    
    %Preallocation 
    u = zeros(3,size(Sn,1));                        %Control law  

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
        e(2) = J-Jref;                            %Jacobi Constant error
        e = e.';
       
        %Control law  
        u(:,i) = -real(K*e);           
    end
end