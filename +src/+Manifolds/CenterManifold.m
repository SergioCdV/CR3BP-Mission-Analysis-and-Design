%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/01/25
% File: CenterManifold.m 
% Issue: 0 
% Validated: 

%% Center Invariant Manifold %%
% This class definition provides the implementation of a general center
% invariant manifold

classdef CenterManifold < src.InvariantManifold
    properties
        Branch;
        eps;
    end

    methods
        % Basic constructor
        function [obj] = CenterManifold()
            % Parent constructor 
            obj@src.InvariantManifold();
        end
    end
end