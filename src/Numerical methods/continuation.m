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

% New versions: 

function [x, state] = continuation(x0, num, method, parameter, parameter_value, object, setup, varargin)
    %Select method to continuate the initial solution
    switch (method)
        case 'SPC' 
            [x, state] = SP_continuation(x0, num, parameter, parameter_value, object, setup, varargin);
        case 'PAC' 
            [x, state] = PA_continuation(x0, num, parameter, parameter_value, object, setup, varargin);
        otherwise 
            disp('No valid option was selected.');
            x = [];
            state = false;
    end
end

%% Auxiliary function 
function [x, state] = SP_continuation(x0, num, parameter, parameter_value, object, setup, varargin)
    switch (object)
        case 'Orbit'
            %Constants 
            state_dim = 6;                          %Phase space dimension 
            bifValue = 1;                           %Value of the Henon stability index to detect a bifurcation
            mu = setup(1);                          %Reduce gravitational parameter of the system
            n = setup(2);                           %Maximum number of iteratinos for the differential correction process 
            tol = setup(3);                         %Tolerance for the differential correction method
            ds = 1e-8;                              %Continuation step (will vary depending on the solution stability)   
            GoOn = true;                            %Boolean to stop the continuation process
            i = 1;                                  %Continuation iteration
            Phi = eye(state_dim);                   %Initial STM 
            Phi = reshape(Phi, [1 state_dim^2]);    %Reshape initial STM 
            
            %Preallocate solution
            step = [ds zeros(1,num-1)];             %Family continuation vector
            state = zeros(1,num);                   %Preallocate convergence solution
            x = zeros(num, state_dim+state_dim^2);  %Preallocate initial seeds
            T = zeros(1,num);                       %Period of each orbit
            s = zeros(1,num);                       %Stability index of each orbit
            x(1,1:state_dim) = x0+step;             %Initial solution
            x(1,state_dim+1:end) = Phi;             %Append initial STM
            
            %Main computation
            while (i <= num) && (GoOn)
                switch (parameter)
                    case 'Energy'                      
                        %Differential correction
                        [Y, state(i)] = differential_correction('Periodic MS', mu, shiftdim(x(i,:)), n, tol, varargin);
                        STM = reshape(Y.Trajectory(state_dim+1:end), state_dim, state_dim); 
                        
                        %Study stability 
                        [s(i), stm_state] = henon_stability(STM); 
                        bif_flag = (abs(bifValue-s(i)) == 0);
                        
                        %Compute the energy of the solution 
                        C = jacobi_constant(mu, shiftdim(Y.Trajectory(end,1:state_dim)));
                        if (isnan(parameter_value))
                            parameter_value = 10^15;
                        else
                            par_error = abs(parameter_value-C);
                        end
                        
                        %Convergence and stability analysis
                        if (stm_state) && (bif_flag) && (par_error > tol)   
                            iter = iter+1;                                          %Update iteration value
                            T(i) = Y.Period;                                        %Update the period vector
                            x(i+1,1:state_dim) = Y.Trajectory(1,1:state_dim)+step;  %Update initial conditions
                            
                            %Modify the step depending on the proximity to a bifurcation
                            step(1) = par_error*abs(bifValue-s(i))*step(1);
                        else
                            GoOn = false;                                           %Stop the process
                        end  
                                           
                    case 'Period'
                end      
            end
            
            %Output 
            Output.Stability = s;        
            Output.Seeds = x(:,1:state_dim);   
            Output.Period = T;   
            
        case 'Torus'
        otherwise
            x = []; 
            state = false;
            disp('No valid object was selected'); 
            disp(' '); 
    end
end

function [x, state] = PA_continuation(x0, num, parameter, parameter_value, object, setup, varargin) 
end