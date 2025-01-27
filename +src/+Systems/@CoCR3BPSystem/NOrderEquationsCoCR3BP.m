%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 28/12/24
% File: NOrderEquationsCoCR3BP.m 
% Issue: 0 
% Validated: 

%% N-th order Equations of the Co-orbital CR3BP Dynamics %%

% Inputs: 
% Outputs: - vector ds, the differential vector field

% New versions: 

function [ds] = NOrderEquationsCoCR3BP(t, j, s, u, params)
    % Define the initial phase space vector
    r_t = s(1:3,:);                                 % Target synodic position vector
    s_r = s(7:12,:);                                % Relative synodic state vector
    rho = s_r(1:3,:);                               % Relative synodic position vector
    rho_norm = dot(rho, rho, 1);                    % Norm of the relative position vector
    D = sqrt( dot(r_t, r_t, 1) );                   % Distance of the reference to the nearest primary

    % Relevant system parameters
    mu = params(1);                                 % Gravitational parameter of the system
    N = params(end);                                % Order of the approximation

    if ( length(params) == 2 )
        % Coefficients of the Hamiltonian expansion
        cn = src.Systems.CoCR3BPSystem.CoLegendreCoefficients(mu, r_t, N);

    else
        % Additional parameters
        L = params(2);                                  % Index of the reference libration point
        gamma = params(3);                              % Distance of the libration point to the closest primary
        
        % Coefficients of the Hamiltonian expansion
        cn = src.Systems.CR3BPSystem.LegendreCoefficients(mu, L, gamma, N);
    end

    % Linear acceleration 
    O = zeros(3);                                   % 3 by 3 null matrix
    I = eye(3);                                     % 3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];                 % Hat map dyadic of the angular velocity for the synodice reference frame
    A = [O I; O 2 * Omega];                         % Constant state matrix 

    ds = A * s_r;
    ds(4:6,:) = ds(4:6,:) + [rho(1:2,:); zeros(1,size(s_r,2))];

    grad = zeros(3, N + 1);                         % Pre-allocation of the gradient vector 
    Tn = zeros(1, N + 1);                           % Pre-allocation of the Legendre polynomials

    for j = 1:size(s,2)
        % Initial values
        prod_dot = dot(r_t(:,j), rho(:,j), 1);
        
        Tn(1,1) = 1 / D(j); 
        Tn(1,2) = prod_dot / D(j)^3;

        grad(:,1) = [0 0 0].'; 
        grad(:,2) = r_t(:,j).' / D(j)^3; 

        for i = 3:(N+1)
            n = i - 1;      % Order of the polynomial

            % New Legendre polynomials
            Tn(i) = (2 * n - 1) / n * prod_dot * Tn(1,i-1) - (n - 1) / n * rho_norm(j) * Tn(1,i-2);
            Tn(i) = Tn(i) / D(j);

            % New gradients 
            grad(:,i) = (2 * n - 1) / n * ( r_t(:,j) * Tn(1,i-1) + prod_dot * grad(:,i-1) );
            grad(:,i) = grad(:,i) - (n - 1) / n * ( rho_norm(j) * grad(:,i-2) + 2 * Tn(1,i-2) * rho(:,j) );

            grad(:,i) = grad(:,i) / D(j);
        end

        acc = sum(cn .* grad, 2);
        ds(4:6,j) = ds(4:6,j) + acc; 
    end
 
    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;
end