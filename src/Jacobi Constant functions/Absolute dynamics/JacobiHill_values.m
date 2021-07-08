%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 30/09/20
% File: JacobiHill_values.m 
% Issue: 0 
% Validated: 

%% Jacobi-Hill value %%
% For a given gravitational parameter mu, and a libration point position vector, this 
% function computes the critical values that opens the four Hill's regions (the 
% four energy cases/realms for which motion is possible).

% Inputs: - scalar mu, the reduced central bodies gravitational parameter.

% Outputs: - the scalar J, the associated Jacobi constants that opens the Hill region.

% New versions: 

function [J] = JacobiHill_values(mu)
    %Null velocity vector (definition of Jacobi-Hill critical value) 
    v = [0; 0; 0]; 
    
    %Obtain libration points coordinates 
    L = libration_points(mu);
    
    %Preallocation 
    J = zeros(size(L,2),1);
    
    %Compute critical value 
    for i = 1:size(L,2)
        r = shiftdim(L(1:3,i));         %Libration point position vector
        s = [r; v];                     %Phase state vector
        J(i) = jacobi_constant(mu, s);  %Associated Jacobi Constant
    end
end