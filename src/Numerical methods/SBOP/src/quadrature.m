%% Project: 
% Date: 08/11/22

%% Quadrature scheme %%
% Function to compute the quadrature scheme weights and associated
% integration domain

% Inputs: - structure system, containing the physical information of the
%           2BP of interest
%         - vector initial_coe, the initial orbital elements 
%         - vector final_coe, the final orbital elements 
%         - scalar K, an initial desired revolutions value 
%         - scalar T, the maximum allowed acceleration 
%         - structure setup, containing the setup of the algorithm in general

% Outputs: - vector tau, the domain of integration
%          - vector W, the quadrature weights 
%          - scalar J, the Jacobian of the domain transformation


function [tau, W, J, f, D] = quadrature(n, m, sampling_distribution)
    % Final sampling distribution setup
    switch (sampling_distribution)
        case 'Chebyshev'
            % Sanity check on the quadrature number of points 
            if (m <= max(n))
                m = max(n)+1;
            end

            if (mod(m,2) ~= 0)
                m = m+1;
            end

            % Jacobian domain transformation 
            J = 0.5;
            
            % Clenshaw-Curtis Quadrature weights
            [tau, W] = CC_weights(m);

            % Domain transformation 
            f = @(t0, tf, tau)[(tf+t0)/2+(tf-t0)/2*tau; (tf-t0) * J * ones(1,length(tau))];

        case 'Legendre'
            % Sanity check on the quadrature number of points 
            if ((m-1)/2 <= max(n))
                m = 2*max(n)+1;
            end

            % Jacobian domain transformation 
            J = 0.5;

            % Guass-Legendre Quadrature weights
            [tau, W, D] = LGL_weights(m);

            % Domain transformation 
            f = @(t0, tf, tau)[(tf+t0)/2+(tf-t0)/2*tau; (tf-t0) * J * ones(1,length(tau))];

        case 'Linear'
            % Jacobian domain transformation 
            J = 1;

            % Integration domain 
            tau = sampling_grid(m, sampling_distribution, '');

            % Newton-Cotes Quadrature weights
            W = NC_weights(2);
            W = [];

            % Domain transformation 
            f = @(t0, tf, tau)[(tf-t0) * tau; (tf-t0) * ones(1,length(tau))];

        case 'Normal'
            % Jacobian domain transformation 
            J = 1;

            % Integration domain 
            tau = sampling_grid(m, sampling_distribution, '');

            % Newton-Cotes Quadrature weights
            W = NC_weights(3);
            W = [];

            % Domain transformation 
            f = @(t0, tf, tau)[(tf-t0) * tau; (tf-t0) * ones(1,length(tau))];

        case 'Random'
            % Jacobian domain transformation 
            J = 1;

            % Integration domain 
            tau = sampling_grid(m, sampling_distribution, '');

            % Newton-Cotes Quadrature weights
            W = NC_weights(10);
            W = [];

            % Domain transformation 
            f = @(t0, tf, tau)[(tf-t0) * tau; (tf-t0) * ones(1,length(tau))];

        case 'Newton-Cotes'
            % Jacobian domain transformation 
            J = 1;

            % Integration domain 
            tau = sampling_grid(m, sampling_distribution, '');

            % Newton-Cotes Quadrature weights
            W = NC_weights(10);
            W = [];

            % Domain transformation 
            f = @(t0, tf, tau)[(tf-t0) * tau; (tf-t0) * ones(1,length(tau))];

        case 'Trapezoidal'
            % Jacobian domain transformation 
            J = 1;

            % Integration domain 
            tau = sampling_grid(m, 'Linear', '');

            % No weights
            W = [];

            % Domain transformation 
            f = @(t0, tf, tau)[(tf-t0) * tau; (tf-t0) * ones(1,length(tau))];

        case 'Bernstein'
            % Sanity check on the quadrature number of points 
            if (m <= max(n))
                m = max(n)+1;
            end

            % Jacobian domain transformation 
            J = 1; 

            % Scaled Gauss-Legendre Quadrature weights
            [tau, W, D] = LGL_weights(m);
            W = W/2;
            tau = 0.5*(tau+1);

            % Domain transformation 
            f = @(t0, tf, tau)[(tf-t0) * tau; (tf-t0) * 0.5 * ones(1,length(tau))];
            
        otherwise
            error('No valid quadrature was selected');
    end
end