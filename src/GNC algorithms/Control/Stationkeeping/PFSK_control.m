%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 02/12/21
% File: OPFSK_control.m 
% Issue: 0 
% Validated: 02/12/21

%% Online Primer Floquet Stationkeeping %%
% This script contains the function to compute the control law by means of an PFSK controller.

% Inputs: - vector Sn, the chaser spacecraft state
%         - vector t, the elapsed time
%         - array J, the Floquet exponents matrix
%         - vector lambda, the initial conditions for the primer vector
%         - string cost_function, to optimize on the l1 or l2 norm of the
%           control vector
%         - scalar Tmax, the maximum available thrust for L1 minimization

% Output: - array u, containing the required control law

% New versions: 

function [u] = PFSK_control(t, Sn, P, J, lambda, cost_function, Tmax)
    %Constants 
    m = 6;       %Phase space dimension

    %Preallocation of the control vector 
    u = zeros(3,size(Sn,1));

    %Control input matrix in the synodic frame
    B = [zeros(3,3); eye(3)];

    for i = 1:size(Sn,1)
        %Reassemble the Floquet projection matrix 
        Phi = reshape(Sn(i,m+1:end), [m m]);
        F = Phi*P*expm(-J*t(i));

        %Compute the projected control matrix in the Floquet space 
        V = F^(-1)*B;

        %Switch depending on the cot function to minimize
        switch (cost_function)
            case 'L1' 
                %Primer vector 
                p = -pinv(V)*expm(-J*t(i))*lambda; 

                %Direction of the control vector
                if (norm(p) ~= 0)
                    uv = p/norm(p);
                else
                    uv = zeros(3,1);
                end

                %Final control vector
                delta = 1e6;
                u(:,i) = Tmax/(1+exp(-delta*(norm(p)-1)))*uv;

            case 'L2'
                u(:,i) = -V.'*expm(-J*t(i))*lambda;          

            otherwise
                error('No valid vector norm to be minimized was selected');
        end

        u(:,i) = real(u(:,i));
    end
end