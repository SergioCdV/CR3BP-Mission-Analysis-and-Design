%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: setup_path.m 
% Issue: 0 
% Validated: 

%% Set up path %%
% This scripts provides the function to generate the needed search paths 

function setup_path()
    %Generate the search paths of the main subfolders of the AOMAT program
    mainFolder = fileparts(which('setup_path'));
    folderDyn = strcat(mainFolder,'\Dynamics');
    folderNum = strcat(mainFolder,'\Numerical methods');
    folderInv = strcat(mainFolder,'\Invariant objects');
    folderJaC = strcat(mainFolder,'\Jacobi Constant functions');
    folderRef = strcat(mainFolder,'\Reference frames');
    folderTes = strcat(mainFolder,'\Test');
    
    %Add paths
    addpath(genpath(folderDyn), '-end');
    addpath(genpath(folderNum), '-end');
    addpath(genpath(folderInv), '-end');
    addpath(genpath(folderJaC), '-end');
    addpath(genpath(folderRef), '-end');
    addpath(genpath(folderTes), '-end');
end