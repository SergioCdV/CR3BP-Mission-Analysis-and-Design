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
%         - vector S, the time evolution of the spacecraft state 
%         - scalar t, the time at which the GNC scheme is evaluated
%         - cell array varargin, used whenever a target spacecraft uses GNC

% Output: - guidance law Sg, the computed guidance law
%         - estimated state Sn, the output of the navigation filter
%         - control u, the control law to steer the state to the guidance
%           law requirements

function [Sg, Sn, u] = GNCc_handler(GNC, St, S, t)
    %Extract the algorithms flags 
    guidance = GNC.Algorithms.Guidance;              %Selected the guidance algorithm 
    navigation = GNC.Algorithms.Navigation;          %Selected the navigation algorithm 
    control = GNC.Algorithms.Control;                %Selected the control algorithm

    %Navigation module 
    noise = GNC.Navigation.NoiseVariance; 
    Sn = S;
    if (any(noise ~= 0))
        Sn(1,1:6) = S(1,1:6)+normrnd(0,noise,1,6);
    end
    
    %Guidance module 
    switch (guidance)
        case 'APF'
            %Guidance parameters 
            safe_corridor = GNC.Guidance.APF.Safety;        %Safety parameters
            Penalties = GNC.Guidance.APF.Penalties;         %Guidance core parameters
            So = GNC.Guidance.APF.Obstacles;                %Relative position of the obstacles
            
            %Preallocation 
            Sg = zeros(size(S,1),GNC.Guidance.Dimension);
            
            %Guidance law
            for i = 1:size(S,1)
                Sg(i,:) = APF_guidance(safe_corridor, Penalties, So, t, S(i,:).', false);
            end
            
        case 'CTR'
            %Guidance regression
            order = GNC.Guidance.CTR.Order;                 %Order of the approximation
            TOF = GNC.Guidance.CTR.TOF;                     %Time of flight
            Cp = GNC.Guidance.CTR.PositionCoefficients;     %Coefficients of the Chebyshev approximation
            Cv = GNC.Guidance.CTR.VelocityCoefficients;     %Coefficients of the Chebyshev approximation
            
            u = (2*t-TOF)/TOF;                              %Normalized domain
            T = CH_basis('first', order, u);                %Polynomial basis 
            p = Cp*T;                                       %Position trajectory
            v = Cv*T;                                       %Velocity trajectory

            switch (control) 
                case 'LQR'
                    Ci = GNC.Guidance.CTR.IntegralCoefficients;     %Coefficients of the Chebyshev approximation
                    g = Ci*T;                                       %Integral of the position trajectory
                case 'SDRE'
                    Ci = GNC.Guidance.CTR.IntegralCoefficients;     %Coefficients of the Chebyshev approximation
                    g = Ci*T;                                       %Integral of the position trajectory
                case 'HDRE'
                    Ci = GNC.Guidance.CTR.IntegralCoefficients;     %Coefficients of the Chebyshev approximation
                    g = Ci*T;                                       %Integral of the position trajectory
                otherwise
                    Cg = GNC.Guidance.CTR.AccelerationCoefficients; %Coefficients of the Chebyshev approximation
                    g = Cg*T;                                       %Acceleration trajectory
            end
            
            %Guidance trajectory
            Sg = [p.' v.' g.'];                             
            
        otherwise
            Sg = zeros(size(S,1), GNC.Guidance.Dimension);          %No guidance requirements
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
            s0 = [St(1,1:6).'; S(1,1:6).'];                 %Initial conditions
            
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
            s0 = [St(1,1:6).'; S(1,1:6).'];                 %Initial conditions
            
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
            s0 = [St(1,1:6).'; S(1,1:6).'];                 %Initial conditions
            
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
            u = LQR_control(model, mu, Sg, Sn(:,1:9), target, Ln, gamma, Q, M);
            
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
            u = SDRE_control(model, mu, Sg, Sn, St(:,1:6), Ln, gamma, Q, M);

        case 'HDRE'
            %System characteristics 
            mu = GNC.System.mu;                 %Systems's reduced gravitational parameter
            Ln = GNC.System.Libration(1);       %Libration point ID
            gamma = GNC.System.Libration(2);    %Libration point distance to the neareast primary
            
            %Controller parameters
            model = GNC.Control.H.Model;        %Linear model to use
            SDRE = GNC.Control.H.SDRE_flag;     %Time-varying dynamics flag
            W = GNC.Control.H.Variance;         %Variance of the noise vector
            if (SDRE)
                target = St;                    %Final desired target state
            else
                target = GNC.Control.H.Target;  %Final desired target state
            end
            
            %Control law
            u = HDRE_control(model, mu, Sg, Sn, target, Ln, gamma, SDRE, W);
            
        case 'SMC'
            %System characteristics 
            mu = GNC.System.mu;                             %Systems's reduced gravitational parameter
            
            %Controller parameters 
            parameters = GNC.Control.SMC.Parameters;        %Parameters of the controller
            model = 'Encke';                                %Dynamics vector field approximation model
            Stotal = [St(:,1:6) S(:,1:6)];                  %Complete phase space vector
            
            %Control law
            u = SMC_control(mu, Sg, Stotal, parameters, model);
            
        case 'MPC'
            %System characteristics 
            mu = GNC.System.mu;                            %Systems's reduced gravitational parameter
            
            %Controller parameters 
            TOF = GNC.Control.MPC.TOF;                     %Rendezvous time of flight  
            cost_function = GNC.Control.MPC.Cost;          %Cost function for the differential corrector scheme
            Tmin = GNC.Control.MPC.Thruster(1);            %Thruster error model
            Tmax = GNC.Control.MPC.Thruster(2);            %Thruster error model
            core = GNC.Control.MPC.Core;                   %Optimization core to solve the optimal problem
            method = GNC.Control.MPC.Method;               %Nonlinear method to solve the optimal problem
            
            %Control law
            [~, u, ~] = MPC_control(mu, cost_function, Tmin, Tmax, TOF, St, core, method);
            
        case 'APF'
            %System characteristics 
            mu = GNC.System.mu;                               %Systems's reduced gravitational parameter
            
            %Controller parameters 
            TOF = GNC.Control.MPC.TOF;                        %Rendezvous time of flight  
            safe_corridor = GNC.Control.APF.safe_corridor;    %Safe corridor parameters
            Penalties = GNC.Control.APF.Penalties;            %Controller parameters
            So = GNC.Control.APF.Obstacles;                   %Obstacles position
            
            %Control law
            [~, ~, u] = APF_control(mu, safe_corridor, Penalties, So, TOF, [Sn St]);


        %Relative stationkeeping
        case 'MFSK'
            %Stationkeeping parameters
            mu = GNC.System.mu;                        %Systems's reduced gravitational parameter
            T = GNC.Control.MFSK.Period;               %Period of the target orbit 
            tol = GNC.Control.MFSK.Tolerance;          %Tolerance for the differential corrector process
            constraint = GNC.Control.MFSK.Constraint;  %Constraint boolean for energy tracking
            
            %Stationkeeping control law
            u = MFSK_control(mu, T, Sn, tol, constraint);  

        case 'RFSK'
            %Stationkeeping parameters
            method = GNC.Control.RFSK.method;           %Solver to be used
            mu = GNC.System.mu;                         %Systems's reduced gravitational parameter
            K = GNC.Control.RFSK.K;                     %State error controller
            Q = GNC.Control.RFSK.Q;                     %State error penalty
            R = GNC.Control.RFSK.M;                     %State error penalty
            Hg = GNC.Control.RFSK.Reference;            %Reference energy state
            L = GNC.Control.RFSK.FloquetModes;          %Floquet modes of the reference trajectory
            F = GNC.Control.RFSK.FloquetDirections;     %Floquet modes of the reference trajectory
            
            % Stationkeeping control law
            u = RFSK_control(method, mu, t, L, F, St, Hg, Sn, Q, R, K); 
         
        case 'PFSK'
            %Stationkeeping parameters
            J = GNC.Control.PFSK.FloquetExponents;             %Floquet exponents of the reference trajectory
            lambda = GNC.Control.PFSK.InitialPrimer;           %Floquet modes of the reference trajectory
            P = GNC.Control.PFSK.FloquetDirections;            %Floquet modes of the reference trajectory
            cost_function = GNC.Control.PFSK.CostFunction;     %Cost function to minimize
            T = GNC.Control.PFSK.Period;                       %Orbital period of the target orbit

            switch (cost_function)
                case 'L1'
                    Tmax = GNC.Control.PFSK.MaxThrust;         %Maximum available thrust
                otherwise
                    Tmax = 0;                                  %Maximum available thrust
            end

            %Stationkeeping control law
            u = PFSK_control(t, T, Sn, P, J, lambda, cost_function, Tmax); 

        case 'LSB'
           % System characteristics 
           mu = GNC.System.mu;                 % Systems's reduced gravitational parameter
           Ln = GNC.System.Libration(1);       % Libration point ID
           gamma = GNC.System.Libration(2);    % Libration point distance to the neareast primary
           
           % Control law
           u = zeros(3,size(Sn,1));
           for i = 1:size(Sn,1)
                % Time of flight
                switch (GNC.LSB.Method) 
                   case 'Dynamics shape-based'
                       tf = 1e-2;
                   otherwise
                       tf = GNC.LSB.Parameters.TOF-t(i);
                end

                if (tf ~= 0)
                    [~, u_c, ~, ~] = LSB_guidance(mu, Ln, gamma, Sn(i,:), GNC.LSB.Method, tf, GNC.LSB.Parameters); 
                    u(:,i) = u_c(:,1);  
                else
                    a = 1;
                end
           end

        otherwise
            warning('No valid control law was selected')
            u = zeros(GNC.Control.Dimension,size(Sn,1));    %No control requirements
    end
end