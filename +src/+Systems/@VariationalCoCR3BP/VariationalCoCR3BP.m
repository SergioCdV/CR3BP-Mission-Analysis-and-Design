%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 28/12/24
% File: VariationalCoCR3BP.m 
% Issue: 0 
% Validated:

%% Variational Equations of the co-orbital CR3BP %% 
% This class implements the definition of the dynamics of the variational
% system of equations for the co-orbital CR3BP

classdef VariationalCoCR3BP < src.DynamicalSystems.ContinuousSystem

    methods
        % Constructor of the class
        function [obj] = VariationalCoCR3BP(n)
           % Constructor of the super class
           obj@src.DynamicalSystems.ContinuousSystem(n^2, 0);
           
           % Dynamics of the problem
           obj.Dynamics = @(t, j, s, u, params)src.Systems.VariationalCoCR3BP.VariationalEquationsCoCR3BP(t, j, s, params);
        end
    end

    methods (Static)
        [ds] = VariationalEquationsCoCR3BP(t, j, s, params);        % First order variational equations of the CR3BP
    end
end