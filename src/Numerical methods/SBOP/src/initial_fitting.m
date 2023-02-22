%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/02/22

%% Initial fitting %%
% Function to estimate the trajectory approximation

% Inputs: - scalar n, the degree of the approximation 
%         - vector tau, the control parameter vector 
%         - C, the initial trajectory estimation
%         - string basis, the Bernstein polynomial basis to be used

% Outputs: - array P, the estimation of the boundary control points as a
%            cell
%          - array C, the initial estimation of the spacecraft state vector

function [P, C] = initial_fitting(Problem, basis, tau, s)
    % Constants 
    n = Problem.PolOrder;       % Polynomial order for each state variable 
    L = Problem.DerDeg;         % Order of the dynamics (maximum derivative order)
    
    % Preallocation of the control points and the expansion polynomials
    P = zeros(length(n), max(n)+1); 
    B = state_basis(n, L, basis, tau);

    % Compute the position control points leveraging the complete state vector
    if (L > 1)
        C = [s(1:Problem.StateDim,:) s(1+Problem.StateDim:2*Problem.StateDim,:)];
    else
        C = s(1:Problem.StateDim,:);       
    end

    for i = 1:length(n)
        if (L > 1)
            A = [B{i}(1:n(i)+1,:) B{i}(n(i)+2:2*(n(i)+1),:)];
        else
            A = B{i}(1:n(i)+1,:);
        end
        P(i,1:n(i)+1) = C(i,:)*pinv(A);
    end
    
    % Evaluate the state vector
    C = evaluate_state(P, B, n, L);
end