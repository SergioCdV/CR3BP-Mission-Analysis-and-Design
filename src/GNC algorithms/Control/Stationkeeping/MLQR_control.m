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

function [u] = MLQR_control(mu, t, T, L, F, Sg, Sn, Q, M)
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

    switch (model)
        case 'Floquet'
            %Compute the state matrix
            A = L;
        
            for i = 1:size(Sn,1)
                %Compute the actual state
                Monodromy = reshape(Sn(i,n+1:n+n^2), [n n]);
        
                E = Monodromy*F*expm(-L*t(i));      %Floquet basis                                                                                                           
                invE = E^(-1);                      %Floquet change of basis
                X = invE*dSn(i,1:n).';              %Floquet constant coefficients
                e = X-Sg(i,n+4:end).';              %State error
        
                %Compute the gradient of the integrals 
                V = invE*B;
                [K,~,~] = lqr(A,V,Q,M);
                u(:,i) = real(-K*e);
            end
        otherwise 
            error('No valid linear model was selected')
    end
end