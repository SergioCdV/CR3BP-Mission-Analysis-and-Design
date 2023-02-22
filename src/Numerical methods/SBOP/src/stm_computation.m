%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% State Transition Matrix Computation %% 
% Function to compute the STM of the problem

% Inputs: - scalar mu, the gravitational parameter of the central body  
%         - array St, defining the target trajectory in time
%         - array s, the optimal relative trajectory
%         - vector tau, the vector of collocation points
%         - string method, the method to compute the STM

% Outputs: - array STM, the evolution of the STM in time  

function [STM] = stm_computation(mu, St, s, tau, method)
    switch (method)
        case 'Analytical'
            [STM] = analytical_STM(n, s, basis);
            
        case 'Numerical'
             % Integration
            options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);     % Integration setup
            Phi0 = reshape(eye(6), [36 1]);                            % Initial conditions
    
            [~, STM] = ode113(@(t,s)var_equations(mu, tf, St, n, P, sampling_distribution, basis, t, s), tau, Phi0, options);
            STM = STM.'; 
    end
end

%% Auxiliary functions 
% Function to analytically compute the STM 
function [STM] = analytical_STM(n, s, basis)
    % Constants 
    m = 6;          % Phase space dimension

    % Initialization
    P = zeros(m^2/2,size(P0,2));  

    % Switch the polynomial basis to be used
    switch (basis)
        case 'Bernstein'                
            % Control points for a nonorthogonal Bézier curve
            P(:,1) = reshape([eye(3) zeros(3)], [m^2/2 1]);
            P(:,2) = reshape([eye(3) diag(1./n)], [m^2/2 1]);

        otherwise
            P = P0; 
            % Compute the partial state evolution 
            N = size(P,1);                   % Number of state variables
            C = zeros(2*N,2);                % Preallocation for speed
            for i = 1:length(n)
                C(i,:) = P0(i,3:n(i)-1)*B{i}(3:n(i)-1,[1 end]);
                C(length(n)+i,:) = P0(i,3:n(i)-1)*B{i}(n(i)+4:2*n(i),[1 end]);
            end

            % Assemble the linear system 
            C = [x0.'-C(:,1) xf.'-C(:,end)];
            C = [C(1:size(P,1),:) C(size(P,1)+1:2*size(P,1),:)];
            
            % Control points for an orthogonal Bézier curve
            for i = 1:length(n)
                index = [1 2 n(i) n(i)+1 n(i)+2 n(i)+3 2*n(i)+1 2*(n(i)+1)];
                A = B{i}(index,[1,end]);
                A = [A(1:4,:) A(5:8,:)];
                P(i,[1 2 n(i) n(i)+1]) = C(i,:)*A^(-1);
            end
    end

    % Compute the STM 
    k = 1; 
    STM = zeros(m^2, size(B{1},2));
    for i = 1:size(P,1)
        l = n(k)+1;
        for j = 1:2
            STM(i+size(P,1)*(j-1),:) = P(i,1:l)*B{k}(1+l*(j-1):l*j,:);
        end
        if (mod(i,m) == 0)
            k = k+1; 
        end
    end
end