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
    guidance = GNC.Algorithms.Guidance;              %Selected the guidance algorithm 
    navigation = GNC.Algorithms.Navigation;          %Selected the navigation algorithm 
    control = GNC.Algorithms.Control;                %Selected the control algorithm
    
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
    
    %Control module
    switch (control)    
        %Stationkeeping module 
        case 'HSK'
            %Stationkeeping parameters
            Q = GNC.Control.HSK.Q;                   %State error penalty
            R = GNC.Control.HSK.M;                   %State error penalty
            mu = GNC.System.mu;                      %Systems's reduced gravitational parameter
            u = HSK_control(mu, Sg, Sn, Q, R);       %Stationkeeping control law

        case 'MLQR'
            %Stationkeeping parameters
            Q = GNC.Control.MLQR.Q;                  %State error penalty
            R = GNC.Control.MLQR.M;                  %State error penalty
            mu = GNC.System.mu;                      %Systems's reduced gravitational parameter
            Hg = GNC.Control.MLQR.Reference;         %Reference energy state
            Sg = [Sg repmat(Hg.', size(Sg,1), 1)];   %Reference orbital and energy state
            T = GNC.Control.MLQR.Period;             %Period of the target orbit 
            L = GNC.Control.MLQR.FloquetModes;       %Floquet modes of the reference trajectory
            F = GNC.Control.MLQR.FloquetDirections;  %Floquet modes of the reference trajectory
            
            %Stationkeeping control law
            u = MLQR_control(mu, t, T, L, F, Sg, Sn, Q, R);   

         case 'MFSK'
            %Stationkeeping parameters
            mu = GNC.System.mu;                        %Systems's reduced gravitational parameter
            Jref = GNC.Control.MFSK.Reference;         %Reference energy state
            T = GNC.Control.MFSK.Period;               %Period of the target orbit 
            tol = GNC.Control.MFSK.Tolerance;          %Tolerance for the differential corrector process
            constraint = GNC.Control.MFSK.Constraint;  %Constraint boolean for energy tracking
            
            %Stationkeeping control law
            u = MFSK_control(mu, T, Sn, tol, constraint, Jref);  

        %Artificial objects control law
        case 'TAHO'
            %System characteristics 
            mu = GNC.System.mu;                      %Systems's reduced gravitational parameter
            Sc = Sn;                                 %Artificial halo orbit equilibrium position
            
            %Control law
            u = TAHO_control(mu, Sc);
            
        otherwise 
            u = zeros(6,1);                          %Stationkeeping control law
    end
end