%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: LinearEquationsCoCR3BP.m 
% Issue: 0 
% Validated: 

%% Linear Equations of the Co-orbital CR3BP Dynamics %%
% This function contains the description of the co-orbital CR3BP dynamics vector field

% Inputs: 
% Outputs: - vector ds, the differential vector field

% New versions: 

function [ds] = LinearEquationsCoCR3BP(t, j, s, u, params)
    % Define the initial phase space vector
    r_t = s(1:3,:);                                 % Target synodic position vector
    s_r = s(7:12,:);                                % Relative synodic state vector
    rho = s(7:9,:);                                 % Relative synodic position vector
    
    % Relevant system parameters
    mu = params(1);                                                         % Gravitational parameter of the system
    mup(1) = 1 - mu;                                                        % First primary normalized position
    mup(2) = mu;                                                            % Second primary normalized position
    R(:,1) = reshape(params(2:4), [], 1);                                   % Position vector of the first primary
    R(:,2) = reshape(params(5:7), [], 1);                                   % Position vector of the second primary
    order_flag = params(8);                                                 % Flag to include second order terms for recursive purposes
    
    % Relative position between the primaries and the target 
    Rr(1:3,:) = R(:,1) - r_t;                                               % Position of the target with respect to the first primary
    Rr(4:6,:) = R(:,2) - r_t;                                               % Position of the target with respect to the first primary
    Rr_norm(1,:) = sqrt( dot(Rr(1:3,:), Rr(1:3,:), 1) );
    Rr_norm(2,:) = sqrt( dot(Rr(4:6,:), Rr(4:6,:), 1) );
    ur(1:3,:) = Rr(1:3,:) ./ Rr_norm(1,:);                                  % Unit vector of the relative position of the target with respect to the first primary
    ur(4:6,:) = Rr(4:6,:) ./ Rr_norm(2,:);                                  % Unit vector of the relative position of the target with respect to the second primary
    
    % Relative position unit vector
    rho_norm = sqrt( dot(rho, rho, 1) );
    zero_idx = rho_norm ~= 0;
    u_rho = rho;

    if any(zero_idx)
        u_rho(:,zero_idx) = rho(:,zero_idx) ./ rho_norm(zero_idx);    
    end

    % Relative acceleration vector field
    O = zeros(3);                                   % 3 by 3 null matrix
    I = eye(3);                                     % 3 by 3 identity matrix
    Omega = [0 1 0; -1 0 0; 0 0 0];                 % Hat map dyadic of the angular velocity for the synodice reference frame
    A = [O I; O 2 * Omega];                         % Constant state matrix 

    % Relative acceleration (linear order term) 
    ds = A * s_r;

    for i = 1:length(mup)
        idx = 1 + 3 * (i - 1) : 3 * i;

        % Relative Legendre coefficient          
        c2 = mup(i) ./ Rr_norm(i,:).^3;          
        
        % Compute the acceleration
        cos_theta = dot(u_rho, ur(idx,:), 1);

        for j = 1:size(s,2)
            if (zero_idx(j))
                B = (I - order_flag * rho(:,j) * rho(:,j).' / dot(rho(:,j), rho(:,j)) ).';
            else
                B = I;
            end
            ds(4:6,j) = ds(4:6,j) + c2 * ( (3 * order_flag * cos_theta(j)^2 - 1) * I + 3 * (ur(idx,j) * ur(idx,j).') * B) * rho(:,j);
        end
    end
 
    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;
end