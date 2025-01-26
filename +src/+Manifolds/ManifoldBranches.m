%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 07/01/25
% File: ManifoldBranches.m 
% Issue: 0 
% Validated: 

%% Invariant manifold branches %%
% This class definition provides an enumeration for the different branches
% of the hyperbolic invariant manifolds

classdef ManifoldBranches < double
   enumeration
       Left   (+1)
       Right  (-1)
   end
end