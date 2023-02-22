%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer in the CR3BP using a polynomial shape-based approach

% Inputs: - class Problem, the definition of the optimal control problem
%         - array St, a sampled trajectory of the target state
%         - bool STM_flag, the flag to compute the STM

% Outputs: - array C, the final state evolution matrix
%          - scalar cost, the final cost of the transfer 
%          - array u, a 3xm matrix with the control input evolution  
%          - vector t, the time sampling points final distribution

function [C, cost, u, tau] = RSBOPT_wrapper(Problem, St, STM_flag)
    % Target's orbit high-order approximation
    order = 50;                                             % Order of the polynomial
    theta = linspace(0,2*pi,size(St,1));                    % Anomaly evolution 

    [Cp, Cv, ~] = CTR_guidance(order, theta, St);           % Regression of the torus function
    Cs = [Cp; Cv];                                          % Target's orbit position and velocity coordinates

    Problem.Params = [Problem.Params; reshape(Cs, [], 1)];

    % Optimization
    [C, cost, u, ~, ~, tau, ~, ~] = sb_solver(Problem);

    % Final target trajectory 
    St = target_trajectory(tau, Problem.Params(4), Cs);
        
    % Back transaltion of the origin of the synodic frame 
    C = cylindrical2cartesian(C, true);
    C = [St(1:6,:); C];

    % Integrate the STM 
    if (STM_flag)          
        STM = stm_computation(Problem.Params(1), St, C, tau, 'Numerical'); 
    else
        STM = [];
    end

    % Final evolution
    C = [C; STM];
end