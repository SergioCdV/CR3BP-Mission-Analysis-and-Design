
function [Output, state] = SP_Orbit_continuation(object_number, parametrization, Object, corrector, setup)
    %Constants 
    state_dim = 6;                              %Phase space dimension 
    nodes = 15;                                 %Number of nodes to correct the object

    parameter = parametrization{1};             %Parameter to continuate on the initial object
    parameter_value = parametrization{2};       %Desired parameter value 
    seed = Object{2};                           %Initial orbit seed
    object_period = Object{3};                  %Orbit initial period
    mu = setup(1);                              %Reduce gravitational parameter of the system
    n = setup(2);                               %Maximum number of iterations for the differential correction process 
    tol = setup(3);                             %Tolerance for the differential correction method
    direction = setup(4);                       %Direction to continuate in
            
    %Preallocate solution
    X = zeros(object_number, state_dim);        %Preallocate initial seeds
    T = zeros(1,object_number);                 %Period of each orbit
    stability = zeros(1,object_number);         %Stability index of each orbit
    y = seed;                                   %Initial solution
    
    GoOn = true;                                %Boolean to stop the continuation process
    num = 1;                                    %Continuation iteration
                
    %Main computation
    switch (parameter)
         case 'Energy'  
            %Modify initial conditions 
            ds = direction*(5e-3);                      %Continuation step (will vary depending on the solution stability)
            step = [0 0 ds 0 0 0];                      %Family continuation vector
            y(1,1:state_dim) = y(1,1:state_dim)+step;   %Modify initial conditions 

            %Main loop
            while (num <= object_number) && (GoOn)
               %Differential correction
               [Y, state(num)] = differential_correction(corrector, mu, y, n, tol, nodes, object_period);
               STM = reshape(Y.Trajectory(end,state_dim+1:end), state_dim, state_dim); 

               %Study stability 
               [stability(1:state_dim/2,num), stm_state] = henon_stability(STM); 

               %Compute the energy of the solution 
               C = jacobi_constant(mu, shiftdim(Y.Trajectory(end,1:state_dim)));
               if (isnan(parameter_value))
                   par_error = 1;
               else
                   par_error = abs(parameter_value-C);
               end
               
               %Convergence and stability analysis
               if (stm_state) && (par_error > tol)   
                   T(num) = Y.Period;                            %Update the period vector
                   X(num,:) = Y.Trajectory(1,1:state_dim);       %Save initial conditions

                   %Update initial conditions
                   y = Y.Trajectory(:,1:state_dim);
                   y(1,:) = y(1,:)+step;     
                   num = num+1;                                  %Update iteration value
               else
                   num = num+1;                                  %Update object number
                   %Correct the final desired orbit
                   [Y, state(num)] = differential_correction('Jacobi Constant Multiple Shooting', mu, y, n, tol, ...
                                                             nodes, object_period, Cref);
                   X(num,:) = Y.Trajectory(1,1:state_dim);       %Save initial conditions
                   GoOn = false;                                 %Stop the process
               end  
            end
           
         case 'Period'
            %Modify initial conditions 
            ds = direction*(1e-1);                               %Continuation step 
            object_period = object_period+ds;                    %Modify initial conditions 
            corrector = 'Periodic Multiple Shooting';            %Algorithm corrector

            %Main loop
            while (num <= object_number) && (GoOn)
                %Differential correction
                [Y, state(num)] = differential_correction(corrector, mu, y, n, tol, nodes, object_period);
                STM = reshape(Y.Trajectory(end,state_dim+1:end), state_dim, state_dim); 

                %Study stability 
                [stability(1:state_dim/2,num), stm_state] = henon_stability(STM); 

                %Compute the energy of the solution 
                if (isnan(parameter_value))
                    par_error = 1;
                else
                    par_error = abs(parameter_value-object_period);
                end

                %Convergence and stability analysis
                if (stm_state) && (par_error > tol)   
                    T(num) = Y.Period;                           %Update the period vector
                    X(num,:) = Y.Trajectory(1,1:state_dim);      %Save initial conditions

                    %Update initial conditions
                    y = Y.Trajectory(:,1:state_dim);                            
                    object_period = object_period+ds;            %Update orbit period
                    num = num+1;                                 %Update iteration value
                else
                   num = num+1;                                  %Update object number
                   %Correct the final desired orbit
                   [Y, state(num)] = differential_correction('Periodic Multiple Shooting', mu, y, n, tol, ...
                                                             nodes, object_period);
                   X(num,:) = Y.Trajectory(1,1:state_dim);        %Save initial conditions
                   GoOn = false;                                  %Stop the process
                end  
            end
            
       otherwise
           error('No valid continuation parameter was selected');
    end
   
    %Output         
    Output.Seeds = X;   
    Output.Period = T;
    Output.Stability = stability;        
end