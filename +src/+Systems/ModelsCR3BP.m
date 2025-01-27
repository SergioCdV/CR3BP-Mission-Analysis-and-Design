%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 23/12/24
% File: ModelsCR3BP.m 
% Issue: 0 
% Validated: 

%% Models CR3BP %%
% This class definition provides an enumeration for the different manners
% in which to express the CR3BP dynamics

classdef ModelsCR3BP < double
   enumeration
       Newton (0)           % Classical Newton equations
       Encke  (1)           % Synodic motion as co-orbital to the closest equilibrium point
       OrderN (2)           % Synodic motion as co-orbital to the closest equilibrium point
   end
end