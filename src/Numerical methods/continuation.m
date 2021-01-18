%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/10/20
% File: continuation.m 
% Issue: 0 
% Validated: 

%% Simple continuation %%
% This function contains the algorithm to compute continuation methods for 
% any general dynamical object. 

% Inputs: - scalar num, a number of continuated solutions. 
%         - string object, selecting the dynamical solution to continuate.
%         - F, the vector field defining the dynamics, input as a handle
%           function.
%         - vector x0, an initial valid solution from which continuation
%           will proceed.

% Outputs: - array of solutions x, containing the continuated objects and
%            method state.

% Methods: one-parameter continuation as well as pseudo-arc length
% continuation.

% New versions: 

function [x] = continuation(x0, F, num, object)
    %Select object to continuate 
    switch (object)
        case 'Orbit' 
            x = orbit_continuation(x0, F, ds, num);
        case 'Torus' 
            x = tori_continuation(x0, F, ds, num);
        otherwise 
            disp('No valid option was selected.');
            x = [];
    end
end

%% Auxiliary function 
function [x] = orbit_continuation(x0, F, ds, num)
     %Implement the main loop 
     for i = 1:num
         %Compute constraints vectors for the new initial conditions 
         F0 = F(x0); 
         dX = null(F0);
         x0(i+1,:) = x0(i,:)+ds*dX; 
         
         %Apply the differential correction method
     end  
end

function [x] = tori_continuation(x0, F, ds, num) 
end