%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: setup_path.m 
% Issue: 0 
% Validated: 

%% Set up path %%
% This scripts provides the function to generate the needed search paths 

function setup_path()
    %Generate the search paths of the main subfolders of the program
    mainFolder = fileparts(which('setup_path'));

 	%Add paths
    addpath(genpath(mainFolder), '-end');
end