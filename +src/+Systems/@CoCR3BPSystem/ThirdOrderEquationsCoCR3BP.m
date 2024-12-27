%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: ThirdOrderEquationsCoCR3BP.m 
% Issue: 0 
% Validated: 

%% Third order Equations of the Co-orbital CR3BP Dynamics %%
% This function contains the description of the co-orbital CR3BP dynamics vector field

% Inputs: 
% Outputs: - vector ds, the differential vector field

% New versions: 

function [ds] = ThirdOrderEquationsCoCR3BP(t, j, s, u, params)
    % Define the initial phase space vector
    r_t = s(1:3,:);                                    % Target synodic position vector
    rho = s(7:9,:);                                    % Relative synodic position vector
    
    % Relevant system parameters
    mu = params(1);                                    % Gravitational parameter of the system
    mup(1) = 1 - mu;                                   % First primary normalized position
    mup(2) = mu;                                       % Second primary normalized position
    R(:,1) = reshape(params(2:4), [], 1);              % Position vector of the first primary
    R(:,2) = reshape(params(5:7), [], 1);              % Position vector of the second primary
    order_flag = params(8);                            % Flag to include second order terms for recursive purposes

    Rr(1:3,:) = r_t - R(:,1);                          % Synodic relative position of the target to the first primary
    Rr(4:6,:) = r_t - R(:,2);                          % Synodic relative position of the target to the second primary
    rho_norm = sqrt( dot(rho, rho, 1) );               % Relative distance
    zero_idx = rho_norm ~= 0;

    % Relative acceleration (linear and second order term)
    linear_params = [params(1:end-1); 1];
    ds = src.Systems.CoCR3BPSystem.SecondOrderEquationsCoCR3BP(t, j, s, zeros(3,size(s,2)), linear_params);

    % Relative acceleration (third order term)    
    for i = 1:length(mup)
        % Relative Legendre coefficients 
        idx = 1 + 3 * (i-1) : 3 * i;
        Rr_norm = sqrt( dot(Rr(idx,:), Rr(idx,:), 1) );

        c4 = mup(i) ./ Rr_norm.^5;                     % Third order relative Legendre coefficient

        cos_theta = dot(rho, Rr(idx,:), 1);

        if any(zero_idx)
            cos_theta(zero_idx) = cos_theta(zero_idx) ./ (rho_norm(zero_idx) .* Rr_norm(zero_idx));
        end

        for j = 1:size(s,2)    
            % Compute the acceleration
            Sigma = rho_norm(j).^2 * eye(3) - rho(:,j) * rho(:,j).';
            ds(4:6,j) = ds(4:6,j) + (c4/8) * (4 * rho_norm(j)^2 * (35 * order_flag * cos_theta(j)^4 - 30 *cos_theta(j)^2 + 3) * rho ...
                                              + rho_norm(j) * (140 * cos_theta(j)^3 - 60 * cos_theta(j)^2) * Sigma.' * Rr(idx,j) / Rr_norm(j) );
        end
    end
 
    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;
end