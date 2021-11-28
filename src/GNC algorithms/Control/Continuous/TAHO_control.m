%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 01/11/21
% File: TAHO_control.m 
% Issue: 0 
% Validated: 01/11/21

%% Low Thrusted Artificial Halo Orbits Control %%
% This script contains the function to compute the control law by means of an TAHO controller.

% Inputs: - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - array Sn, the system state

% Output: - vector u, the computed control law

% New versions: 

function [u] = TAHO_control(mu, Sn)
    %Preallocation 
    u = zeros(3,size(Sn,1));                            %Control vector
    
    for i = 1:size(Sn,1)
        %Compute the gradient of the Jacobi constant
        dC = jacobi_gradient(mu, Sn(i,:).');

        %Compute the acceleration vector 
        u(:,i) = (1/2)*dC(1:3).';
    end
end