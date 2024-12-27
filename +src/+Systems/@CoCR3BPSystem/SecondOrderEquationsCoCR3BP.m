%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 27/12/24
% File: SecondOrderEquationsCoCR3BP.m 
% Issue: 0 
% Validated: 

%% Second order Equations of the Co-orbital CR3BP Dynamics %%
% This function contains the description of the co-orbital CR3BP dynamics vector field

% Inputs: 
% Outputs: - vector ds, the differential vector field

% New versions: 

function [ds] = SecondOrderEquationsCoCR3BP(t, j, s, u, params)
    % Define the initial phase space vector
    r_t = s(1:3,:);                                 % Target synodic position vector
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

    rho_norm = sqrt( dot(rho, rho, 1) ); 
    zero_idx = rho_norm ~= 0;

    % Relative acceleration (linear order term)
    linear_params = [params(1:end-1); 1];
    ds = src.Systems.CoCR3BPSystem.LinearEquationsCoCR3BP(t, j, s, zeros(3,size(s,2)), linear_params);                       

    % Relative acceleration (third order term)    
    for i = 1:length(mup)
        % Relative Legendre coefficients  
        idx = 1 + 3 * (i-1) : 3 * i;
        Rr_norm = sqrt( dot(Rr(idx,:), Rr(idx,:), 1) );
        
        c3 = mup(i) ./ Rr_norm.^4;      % Third order relative Legendre coefficient
        
        cos_theta = dot(rho, Rr(idx,:), 1);

        if any(zero_idx)
            cos_theta(zero_idx) = cos_theta(zero_idx) ./ (rho_norm(zero_idx) .* Rr_norm(zero_idx));
        end

        for j = 1:size(s,2)
            % Compute the acceleration
            Sigma = dot(rho(:,j), rho(:,j), 1) * eye(3) - rho(:,j) * rho(:,j).';
            ds(4:6,j) = ds(4:6,j) + (c3/2) * ( 3 * rho_norm(j) * (5 * order_flag * cos_theta(j)^3 - 3 * cos_theta(j)) * rho(:,j) + (15 * cos_theta(j)^2 - 3) * Sigma.' * Rr(idx,j) / Rr_norm(j) );
        end
    end
     
    % Control force 
    ds(4:6,:) = ds(4:6,:) + u;
end