%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/10/20
% File: continuation.m 
% Issue: 0 
% Validated: 

%% Simple continuation %%
% This function contains the algorithm to compute continuation methods for 
% any general dynamical object. 

% Inputs: - vector x0, an initial valid solution from which continuation
%           will proceed.
%         - scalar num, a number of continuated solutions to be computed. 
%         - string method, selecting single parameter continuation (SPC) of pseudo-arc continuation (PAC).
%         - string parameter, for SPC, selecting the parameter to continuate on.
%         - scalar parameter_value, the maximum value of the parameter to stop the process. Maybe NaN.
%         - string object, selecting the dynamical solution to continuate.
%         - vector setup, defining the general performance criteria. Depends on method. 
%         - vector varargin, to be used within a multiple shooting differential corrector as seein in differential_correction.

% Outputs: - array of solutions x, containing the seeds for the required objects and their convergence state.

% Methods: single-parameter continuation as well as pseudo-arc length continuation.

% New versions: Period continuation is not working.

function [x, state] = continuation(object_number, method, parametrization, Object, corrector, setup)
    %Select method to continuate the initial solution
    switch (method)
        case 'SPC' 
            [x, state] = SP_continuation(object_number, parametrization, Object, corrector, setup);
        case 'PAC' 
            [x, state] = PA_continuation(object_number, parametrization, Object, setup);
        otherwise 
            disp('No valid option was selected.');
            x = [];
            state = false;
    end
end

%% Auxiliary function 
function [Output, state] = SP_continuation(object_number, parametrization, Object, corrector, setup)
    %Inputs 
    object = Object{1};     %Type of object to continuate
    seed = Object{2};       %Seed to continuate the object
    
    %Main procedure
    switch (object)
        case 'Orbit'
            %Constants 
            state_dim = 6;                              %Phase space dimension 
            nodes = 15;                                 %Number of nodes to correct the object
            bifValue = 1;                               %Value of the Henon stability index to detect a bifurcation
            
            parameter = parametrization{1};             %Parameter to continuate on the initial object
            parameter_value = parametrization{2};       %Desired parameter value 
            object_period = Object{3};                  %Orbit initial period
            mu = setup(1);                              %Reduce gravitational parameter of the system
            n = setup(2);                               %Maximum number of iterations for the differential correction process 
            tol = setup(3);                             %Tolerance for the differential correction method
            direction = setup(4);                       %Direction to continuate for
  
            GoOn = true;                                %Boolean to stop the continuation process
            i = 1;                                      %Continuation iteration
            
            %Preallocate solution
            state = zeros(1,object_number);             %Preallocate convergence solution
            X = zeros(object_number, state_dim);        %Preallocate initial seeds
            T = zeros(1,object_number);                 %Period of each orbit
            stability = zeros(1,object_number);         %Stability index of each orbit
            y = seed;                                   %Initial solution
                
            %Main computation
            switch (parameter)
                case 'Energy'  
                    %Modify initial conditions 
                    ds = direction*1e-3;                      %Continuation step (will vary depending on the solution stability)
                    step = [ds zeros(1,state_dim-1)];         %Family continuation vector
                    y(1,1:state_dim) = y(1,1:state_dim)+step; %Modify initial conditions 
                    
                    %Main loop
                    while (i <= object_number) && (GoOn)
                        %Differential correction
                        [Y, state(i)] = differential_correction(corrector, mu, y, n, tol, nodes, object_period);
                        STM = reshape(Y.Trajectory(end,state_dim+1:end), state_dim, state_dim); 
                        
                        %Study stability 
                        [stability(1:state_dim/2,i), stm_state] = henon_stability(STM); 
                        
                        %Compute the energy of the solution 
                        C = jacobi_constant(mu, shiftdim(Y.Trajectory(end,1:state_dim)));
                        if (isnan(parameter_value))
                            par_error = 1;
                        else
                            par_error = abs(parameter_value-C);
                        end
                        
                        %Convergence and stability analysis
                        if (stm_state) && (par_error > tol)   
                            T(i) = Y.Period;                              %Update the period vector
                            X(i,:) = Y.Trajectory(1,1:state_dim);         %Save initial conditions
                                                        
                            %Update initial conditions
                            y = Y.Trajectory(:,1:state_dim);
                            y(1,:) = y(1,:)+step;     
                            i = i+1;                                      %Update iteration value
                        else
                            GoOn = false;                                 %Stop the process
                        end  
                    end
                     
                    %Output         
                    Output.Seeds = X;   
                    Output.Period = T;
                    Output.Stability = stability;
                                           
                case 'Period'
                    %Modify initial conditions 
                    ds = 1e-3;                        %Continuation step 
                    object_period = object_period+ds; %Modify initial conditions 
                    
                    %Main loop
                    while (i <= object_number) && (GoOn)
                        %Differential correction
                        [Y, state(i)] = differential_correction('Periodic MS', mu, y, n, tol, nodes, object_period);
                        STM = reshape(Y.Trajectory(end,state_dim+1:end), state_dim, state_dim); 
                        
                        %Study stability 
                        [stability(i), stm_state] = henon_stability(STM); 
                        bif_flag = (abs(bifValue-stability(i)) == 0);
                        
                        %Compute the energy of the solution 
                        if (isnan(parameter_value))
                            par_error = 1;
                        else
                            par_error = abs(parameter_value-object_period);
                        end
                        
                        %Convergence and stability analysis
                        if (stm_state) && (~bif_flag) && (par_error > tol)   
                            T(i) = Y.Period;                              %Update the period vector
                            X(i,:) = Y.Trajectory(1,1:state_dim);         %Save initial conditions
                            
                            %Modify the step depending on the proximity to a bifurcation
                            %step(1) = par_error * mod(abs(bifValue-s(i)), bifValue) * step(1);
                            
                            %Update initial conditions
                            y = Y.Trajectory(:,1:state_dim);                            
                            object_period = object_period+ds;             %Update orbit period
                            i = i+1;                                      %Update iteration value
                        else
                            GoOn = false;                                 %Stop the process
                        end  
                    end
                     
                    %Output         
                    Output.Seeds = X;   
                    Output.Period = T;
                    Output.Stability = stability;

            otherwise 
                Output = []; 
                state = false;     
            end
               
        case 'Torus'
            Output = []; 
            state = false;
            
        otherwise
            Output = []; 
            state = false;
            disp('No valid object was selected'); 
            disp(' '); 
    end
end

function [x, state] = PA_continuation(object_number, parametrization, Object, setup) 
end