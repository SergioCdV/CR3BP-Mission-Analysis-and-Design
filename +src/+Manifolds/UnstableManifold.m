%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 26/01/25
% File: UnstableManifold.m 
% Issue: 0 
% Validated: 

%% Invariant Manifold %%
% This class definition provides the implementation of a general unstable
% invariant manifold

classdef UnstableManifold < src.InvariantManifold
    properties
        Branch;
        eps;
    end

    methods
        % Basic constructor
        function [obj] = UnstableManifold()
            % Parent constructor 
            obj@src.InvariantManifold();
        end
    end
end