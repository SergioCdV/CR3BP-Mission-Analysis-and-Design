%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/01/25
% File: StableManifold.m 
% Issue: 0 
% Validated: 

%% Stable Invariant Manifold %%
% This class definition provides the implementation of a general stable
% invariant manifold

classdef StableManifold < src.InvariantManifold
    properties
        Branch;
        eps;
    end

    methods
        % Basic constructor
        function [obj] = StableManifold()
            % Parent constructor 
            obj@src.InvariantManifold();
        end
    end
end