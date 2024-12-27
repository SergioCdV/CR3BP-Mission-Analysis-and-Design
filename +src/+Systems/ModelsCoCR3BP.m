%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: ModelsCoCR3BP.m 
% Issue: 0 
% Validated: 

%% Models Co-orbital CR3BP %%
% This class definition provides an enumeration for the different manners
% in which to express the co-orbital CR3BP dynamics

classdef ModelsCoCR3BP < double
   enumeration
       Newton     (0)       % Newton form of the non-linear model
       Encke      (1)       % Encke form of the non-linear model
       Linear     (2)       % Linear model
       Order2     (3)       % Second order model
       Order3     (4)       % Third order model
       OrderN     (5)       % N-th order model
       Libration  (6)       % Linear model for the co-orbital problem
       Richardson (7)       % Linear model for the co-orbital problem 
   end
end