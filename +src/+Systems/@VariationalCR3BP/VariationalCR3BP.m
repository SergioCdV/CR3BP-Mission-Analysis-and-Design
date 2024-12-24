%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/12/24
% File: VariationalCR3BP.m 
% Issue: 0 
% Validated:

%% Variational Equations of the CR3BP %% 
% This class implements the definition of the dynamics of the variational
% system of equations for the CR3BP

classdef VariationalCR3BP < src.DynamicalSystems.ContinuousSystem

    methods
        % Constructor of the class
        function [obj] = VariationalCR3BP(varargin)
           % Constructor of the super class
           obj@src.DynamicalSystems.ContinuousSystem(36, 0);
           
           % Dynamics of the problem
           obj.Dynamics = @(t, j, s, u, params)src.Systems.VariationalCR3BP.VariationalEquations(t, j, s, params);
        end
    end

    methods (Static)
        [ds] = VariationalEquationsCR3BP(t, j, s, params);        % First order variational equations of the CR3BP
    end
end