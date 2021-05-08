%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/05/21
% File: apf_guidance.m 
% Issue: 0 
% Validated: 08/05/21

%% Guidance, Navigation and Control Handler %% 
% This script provides the function to activate/desactivate the different
% GNC algorithms, to be used within simulations

% Inputs: - structure GNC, containing the needed flags and parameters for
%           each algorithm

% Output: - guidance law Sg, the computed guidance law
%         - estimated state Sn, the output of the navigation filter
%         - control u, the control law to steer the state to the guidance
%           law requirements

function [Sg, Sn, u] = GNC_handler(GNC)
    %Extract the algorithms flags 
    guidance = GNC.Algorithms.Guidance;         %Selected guidance algorithm 
    navigation = GNC.Algorithms.Guidance;       %Selected navigation algorithm 
    control = GNC.Algorithms.Guidance;          %Selected control algorithm 
    
    %Guidance module 
    switch (guidance)
        case 'APF'
            Sg = apf_guidance();
        otherwise
            Sg = zeros(GNC.Guidance.Dimension,1);           %No guidance requirements
    end
    
    %Navigation module      
    switch (navigation)
        otherwise
            Sn = zeros(GNC.Navigation.Dimension,1);         %No navigation requirements
    end
    
    %Control module
    switch (control)
        case 'LQR' 
            u = LQR_control();
        case 'SDRE'
        case 'SMC'
        otherwise
            u = zeros(GNC.Control.Dimension,1);             %No control requirements
    end
end