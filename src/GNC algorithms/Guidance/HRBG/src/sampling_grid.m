%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 16/04/22

%% Sampling grid %%
% Function to compute the sampling grid

% Inputs: - scalar m, the number of sampling points in the grid
%         - string method, the type of grid to be used
%         - string mode, to select if the sampling grid is based on the
%           the intersection of low-order nodes or the nodes of a
%           high-order polynomial

% Outputs: - vector tau, the sampling points to be used 

function [tau] = sampling_grid(m, method, mode)
    % Sanity check on the distribution mode 
    switch (method)
        case 'Linear'
            mode = '';
        case 'Normal'
            mode = '';
        case 'Random'
            mode = '';
        otherwise
    end

    % Sampling grid generation
    switch (mode)
        case 'Intersection'
            tau = zeros(1,sum(2:m));
            for i = 2:m
                tau(1+sum(2:i-1):sum(2:i)) = grid(i, method);
            end
            tau = unique(tau);
            tau = sort(tau);
        otherwise
            tau = grid(m, method);
    end
end

%% Auxiliary functions 
% Sampling grid computation
function [tau] = grid(m, method)
    switch (method)
        case 'Linear'
            tau = linspace(0,1,m);
        case 'Normal'
            tau = normrnd(0,1,1,m-2);
            tau = sort(tau);
            tau = (tau-min(tau))/(max(tau)-min(tau));
            tau = [0 tau 1];
        case 'Random'
            tau = rand(1,m);
            tau = sort(tau);
        case 'Legendre'
            tau = LG_nodes(m);
        case 'Chebyshev'
            tau = CH_nodes(m);
        case 'Laguerre'
            tau = LR_nodes(m,0);
        case 'Hermite'
            % Depricated
            tau = HT_nodes(m);
        case 'Orthogonal Bernstein'
            % Deprecated
            tau = OB_nodes(m);
        otherwise
            error('An appropriate time array distribution must be specified')
    end
    if (size(tau,1) ~= 1)
        tau = tau.';
    end
end