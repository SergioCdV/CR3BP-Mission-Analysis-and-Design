%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/08/21
% File: GNCt_handler.m 
% Issue: 0 
% Validated: 08/08/21

%% Guidance, Navigation and Control Target Handler %% 
% This script provides the function to activate/desactivate the different
% GNC algorithms for the target spacecraft, to be used within simulations

% Inputs: - structure GNC, containing the needed flags and parameters for
%           each algorithm
%         - vector Sn, the time evolution of the spacecraft state 
%         - scalar t, the time at which the GNC scheme is evaluated

% Output: - guidance law Sg, the computed guidance law
%         - estimated state Sn, the output of the navigation filter
%         - control u, the control law to steer the state to the guidance
%           law requirements

function [Sg, Sn, u] = GNCt_handler(GNC, Sn, t)
    %Extract the algorithms flags 
    guidance = GNC.Algorithms.Guidance;              %Selected guidance algorithm 
    navigation = GNC.Algorithms.Navigation;          %Selected navigation algorithm 
    stationkeeping = GNC.Algorithms.Stationkeeping;  %Selected control algorithm
    
    %Navigation module 
    
    %Guidance module 
    switch (guidance)    
        case 'CTR'
            %Guidance regression
            order = GNC.Guidance.CTR.Order;                 %Order of the approximation
            TOF = GNC.Guidance.CTR.TOF;                     %Time of flight
            Cp = GNC.Guidance.CTR.PositionCoefficients;     %Coefficients of the Chebyshev approximation
            Cv = GNC.Guidance.CTR.VelocityCoefficients;     %Coefficients of the Chebyshev approximation
            Cg = GNC.Guidance.CTR.AccelerationCoefficients; %Coefficients of the Chebyshev approximation
            
            T = zeros(order, length(t));                    %Preallocation of the polynomial basis
            u = (2*t-TOF)/TOF;                              %Normalized domain
            for i = 1:length(u)
                T(:,i) = chebyshev('first', order, u(i));   %Polynomial basis 
            end
            p = Cp*T;                                       %Position trajectory
            v = Cv*T;                                       %Velocity trajectory
            g = Cg*T;                                       %Acceleration trajectory
            Sg = [p.' v.' g.'];                             %Guidance trajectory
            
        otherwise
            Sg = zeros(size(Sn,1), GNC.Guidance.Dimension); %No guidance requirements
    end
    
    %Stationkeeping module 
    switch (stationkeeping)       
        case 'Hamiltonian'
            %Stationkeeping parameters
            Q = GNC.Stationkeeping.HSK.Q;            %State error penalty
            M = GNC.Stationkeeping.HSK.M;            %State error penalty
            mu = GNC.System.mu;                      %Systems's reduced gravitational parameter
            u = HSK_control(mu, Sg, St, Q, R);       %Stationkeeping control law
            
        otherwise 
            u = zeros(6,1);                          %Stationkeeping control law
    end
end