%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 07/01/25
% File: HaloBranches.m 
% Issue: 0 
% Validated: 

%% Halo orbit branches %%
% This class definition provides an enumeration for the different branches
% of the halo orbit family

classdef HaloBranches < double
   enumeration
       Northern  (+1)
       Southern  (-1)
   end
end