%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24 
% File: JacobianCoCR3BP.m 
% Issue: 0 
% Validated: 

%% Jacobian of the co-orbital CR3BP problem %% 
% This script provides a function to compute the jacobian of the co-orbital motion dynamics 
% at a given state

function [J] = JacobianCoCR3BP(mu, s)
    % System parameters 
    mup(1) = 1 - mu;                                    % Reduced gravitational parameter of the first primary
    mup(2) = mu;                                        % Reduced gravitational parameter of the second primary
    R(1:3,1) = [-mu; 0; 0];                             % Synodic position of the first primary
    R(1:3,2) = [1 - mu; 0; 0];                          % Synodic position of the second primary
       
    % State variables 
    r_t = s(1:3,:);                                     % Position of the target 
    rho = s(end - 5 : end - 3,:);                       % Relative position
    
    % Relative position to the primaries 
    r(1:3,:) = r_t - R(:,1);                            % Target position to the first primary
    r(4:6,:) = r_t - R(:,2);                            % Target position to the second primary
    Rr_norm(1,:) = sqrt( dot(r(1:3,:), r(1:3,:), 1) );  % Norm of the relative position of the target to the first primary
    Rr_norm(2,:) = sqrt( dot(r(4:6,:), r(4:6,:), 1) );  % Norm of the relative position of the target to the second primary
    rc(1:3,:) = rho + r(1:3,:);                         % Chaser position to the first primary
    rc(4:6,:) = rho + r(4:6,:);                         % Chaser position to the second primary
    
    % Variational equations
    O = zeros(3);                                       % 3 by 3 null matrix
    I = eye(3);                                         % 3 by 3 identity matrix
    Omega = [0 2 0; -2 0 0; 0 0 0];                     % Coriolis dyadic
    C = [1; 1; 0];                                      % Centrifugal force vector

    % Pre-allocation 
    J = zeros(6, 6 * size(s,2));
    H = zeros(9, size(s,2));                            % Preallocation of the hessian of the potential
    
    q =  zeros(2, size(s,2)); 
    f =  zeros(2, size(s,2));
    df = zeros(2, size(s,2));

    for k = 1:length(mup)
        idx = 1 + 3 * (k - 1) : 3 * k;

        % Encke's variables
        q(k,:) = -dot(2 * r(idx,:) + rho, rho, 1) ./ dot(rc(idx,:), rc(idx,:), 1);                      % Encke variable
        f(k,:) = q(k,:) .* (3 * (1 + q(k,:)) + q(k,:).^2) ./ (1 + (1 + q(k,:)).^(3/2));                 % Encke coefficient
            
        % Derivative of the Encke coefficient
        df(k,:) = -( 3 * sqrt(1 + q(k,:)) ./ dot(rc(idx,:), rc(idx,:), 1) ) .* ( 1 - dot(2 * r(idx,:) + rho, rho) ./ dot(rc(idx,:), rc(idx,:), 1) ); 
    end

    for i = 1:3
        for j = i:3   
            % Derivative of the Encke acceleration field
            gamma = 0;

            % Encke acceleration
            for k = 1:length(mup)
                start = 1 + 3 * (k - 1);
                sq_coeff = rc(start + (i - 1),:) .* rc(start + (j - 1),:);

                if (i == j) 
                    gamma = gamma + (mup(k) ./ Rr_norm(k,:).^3) .* ( (1 + f(k,:)) + df(k,:) .* sq_coeff );
                else
                    gamma = gamma + (mup(k) ./ Rr_norm(k,:).^3) .* ( df(k,:) .* sq_coeff );
                end
            end
                    
            % Hessian of the potential function
            if (i == j)
                H(3 * (i - 1) + j,:) = C(i) - gamma;    % Diagonal terms
            else
                H(3 * (i - 1) + j,:) =    0 - gamma;    % Non-diagonal terms
            end
        end
    end

    % Symmetry constraint
    H(4,:) = H(3,:); 
    H(7,:) = H(3,:); 
    H(8,:) = H(6,:);

    % Compute the first variational equations evaluated at the reference
    for i = 1:size(s,2)
        % Jacobian of the system 
        idx = 1 + 6 * (i - 1) : 6 * i;
        H_idx = reshape(H(:,i), 3, 3);
        J(:,idx) = [O I; H_idx Omega];             
    end                
end