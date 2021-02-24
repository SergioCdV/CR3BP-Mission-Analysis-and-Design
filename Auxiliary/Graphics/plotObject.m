%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/02/21
% File: plotObject.m 
% Issue: 0 
% Validated: 

%% Plot orbit %% 
% This script provides a function to plot objects, in 2D and 3D 

% Inputs: - matrix object containing the trajectory in the phase space (n x m or m x n), 
%           with n the number of timeshots and m the phase space dimension (4 or 6)

function plotObject(object)
    %Analyze the dimension of the object 
    dim = size(object);
    
    %Main procedure
    if (any(dim == 20)) || (any(dim == 42))
        %Select the proper format
        if (dim(1) == 4) || (dim(1) == 6)
            object = object.';
        end
        
        %Plot object depending on its dimension
        if (dim(2) == 4)
            plot(object(:,1), object(:,2));
        else
            plot3(object(:,1), object(:,2), object(:,3));
        end
    else
        disp('No valid object was input. Phase space dimension is not correct.'); 
    end
end