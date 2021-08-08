%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/05/21
% File: GNC_handler.m 
% Issue: 0 
% Validated: 08/05/21

%% Guidance, Navigation and Control Handler %% 
% This script provides the function to activate/desactivate the different
% GNC algorithms, to be used within simulations

% Inputs: - structure GNC, containing the needed flags and parameters for
%           each algorithm
%         - vector St, the time evolution of the target spacecraft
%         - vector Sn, the time evolution of the spacecraft state 
%         - scalar t, the time at which the GNC scheme is evaluated
%         - cell array varargin, used whenever a target spacecraft uses GNC

% Output: - guidance law Sg, the computed guidance law
%         - estimated state Sn, the output of the navigation filter
%         - control u, the control law to steer the state to the guidance
%           law requirements

function [Sg, Sn, u] = GNC_handler(GNC, St, S, t, varargin)
    %Branch between a chaser and a target spacecraft 
    if(~isempty(varargin))
        if (varargin{1})
            [Sg, Sn, u] = GNCt_handler(GNC, S, t);      %Target GNC handler
        else
            error('No valid GNC handler was selected');
        end
    else
        [Sg, Sn, u] = GNCc_handler(GNC, St, S, t);      %Chaser GNC handler
    end
end