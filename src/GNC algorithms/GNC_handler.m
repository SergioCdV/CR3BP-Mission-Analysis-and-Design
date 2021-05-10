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

% Output: - guidance law Sg, the computed guidance law
%         - estimated state Sn, the output of the navigation filter
%         - control u, the control law to steer the state to the guidance
%           law requirements

function [Sg, Sn, u] = GNC_handler(GNC, St, Sn)
    %Extract the algorithms flags 
    guidance = GNC.Algorithms.Guidance;         %Selected guidance algorithm 
    navigation = GNC.Algorithms.Navigation;     %Selected navigation algorithm 
    control = GNC.Algorithms.Control;           %Selected control algorithm
    
    %Navigation module     
    
    %Guidance module 
    switch (guidance)
        case 'APF'
            %Guidance flags 
            dynamics = GNC.Guidance.APF.Dynamics;           %Steady or unsteady APF
            safe_corridor = GNC.Guidance.APF.SafeCorridor;  %Safety corridor boolean flag
            
            %Obstacles 
            obstacles_states = GNC.Guidance.APF.Obstacles;  %Obstacles states
            
            %Simulation parameters 
            dt = GNC.Guidance.TimeStep;                     %Simulation timestep
            
            %Guidance law
            Sg = APF_guidance(dynamics, safe_corridor, Sn(:,1:3), obstacles_states, dt);
            Sg(4:6) = zeros(1,3);
            Sg(7:9) = zeros(1,3);
            
        otherwise
            Sg = zeros(size(Sn,1),GNC.Guidance.Dimension);           %No guidance requirements
    end
        
    %Control module
    switch (control)
        case 'TISS'
            %System characteristics 
            mu = GNC.System.mu;                             %Systems's reduced gravitational parameter
            
            %Controller parameters 
            TOF = GNC.Control.TISS.TOF;                     %Rendezvous time of flight 
            tol = GNC.Control.TISS.Tolerance;               %Differential corrector tolerance 
            cost_function = GNC.Control.TISS.Cost;          %Cost function for the differential corrector scheme
            two_impulsive = GNC.Control.TISS.Impulses;      %Number of impulses for the differential corrector scheme
            
            %Initial conditions 
            s0 = [St(1,1:6).'; Sn(1,1:6).'];                %Initial conditions
            
            switch (cost_function)
                case 'Position'
                   G = Sg(1:3);
                case 'Velocity' 
                   G = Sg(4:6); 
                case 'State'
                   G = Sg; 
                otherwise
                    error('No valid cost function was selected');
            end
            
            %Compute the control scheme
            [~, u, ~] = TISS_control(mu, TOF, s0, tol, cost_function, G, two_impulsive); 
           
        case 'MISS'
            %System characteristics 
            mu = GNC.System.mu;                             %Systems's reduced gravitational parameter
            
            %Controller parameters 
            TOF = GNC.Control.MISS.TOF;                     %Rendezvous time of flight 
            tol = GNC.Control.MISS.Tolerance;               %Differential corrector tolerance 
            cost_function = GNC.Control.MISS.Cost;          %Cost function for the differential corrector scheme
            impulses = GNC.Control.MISS.Impulses;           %Impulses definition for the differential corrector scheme
            
            %Initial conditions 
            s0 = [St(1,1:6).'; Sn(1,1:6).'];                %Initial conditions
            
            %Compute the control scheme
            [~, u, ~] = TISS_control(mu, TOF, s0, tol, cost_function, impulses);
            
        case 'TITA'
            %System characteristics 
            mu = GNC.System.mu;                             %Systems's reduced gravitational parameter
            
            %Controller parameters 
            TOF = GNC.Control.TITA.TOF;                     %Rendezvous time of flight 
            tol = GNC.Control.TITA.Tolerance;               %Differential corrector tolerance 
            cost_function = GNC.Control.TITA.Cost;          %Cost function for the differential corrector scheme
            two_impulsive = GNC.Control.TITA.Impulses;      %Number of impulses for the differential corrector scheme
            penalties = GNC.Control.TITA.Penalties;         %Controller penalties policy
            target_points = GNC.Control.TITA.Targets;       %Target points for stationkeeping
            thruster_model = GNC.Control.TITA.Thruster;     %Thruster error model
            
            %Initial conditions 
            s0 = [St(1,1:6).'; Sn(1,1:6).'];                %Initial conditions
            
            %Compute the control scheme
            [~, u, ~] = TITA_control(mu, TOF, s0, tol, cost_function, two_impulsive, penalties, target_points, thruster_model); 
            
        case 'LQR'     
            %System characteristics 
            mu = GNC.System.mu;                 %Systems's reduced gravitational parameter
            Ln = GNC.System.Libration(1);       %Libration point ID
            gamma = GNC.System.Libration(2);    %Libration point distance to the neareast primary
           
            %Controller parameters
            model = GNC.Control.LQR.Model;      %Linear model to use
            Q = GNC.Control.LQR.Q;              %Penalty matrix on the state error
            M = GNC.Control.LQR.M;              %Penalty matrix on the control effort
            
            target = GNC.Control.LQR.Reference; %Reference position of the target spacecraft
            
            %Control law
            u = LQR_control(model, mu, Sg, Sn(:,1:6), target, Ln, gamma, Q, M);
            
        case 'SDRE'
            %System characteristics 
            mu = GNC.System.mu;                 %Systems's reduced gravitational parameter
            Ln = GNC.System.Libration(1);       %Libration point ID
            gamma = GNC.System.Libration(2);    %Libration point distance to the neareast primary
            
            %Controller parameters
            model = GNC.Control.SDRE.Model;     %Linear model to use
            Q = GNC.Control.SDRE.Q;             %Penalty matrix on the state error
            M = GNC.Control.SDRE.M;             %Penalty matrix on the control effort
            
            %Control law
            u = SDRE_control(model, mu, Sg, Sn(:,1:6), St(:,1:6), Ln, gamma, Q, M);
            
        case 'SMC'
            %System characteristics 
            mu = GNC.System.mu;                             %Systems's reduced gravitational parameter
            
            %Controller parameters 
            parameters = GNC.Control.SMC.Parameters;        %Parameters of the controller
            Stotal = [St(:,1:6) Sn(:,1:6)];                 %Complete phase space vector
            
            %Control law
            u = SMC_control(mu, Sg, Stotal, parameters);
            
        otherwise
            u = zeros(GNC.Control.Dimension,size(Sn,1));    %No control requirements
    end
end