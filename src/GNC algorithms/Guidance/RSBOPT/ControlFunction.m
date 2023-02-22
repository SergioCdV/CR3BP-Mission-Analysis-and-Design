%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(params, beta, t0, tf, t, s)
    % Constants 
    mu = params(1); 
    
    % Computation of the target trajectory
    St = target_trajectory(t(1,:), params(4), reshape(params(5:end), 6, []));

    % Constants of the problem
    mu_r(1) = 1-mu;                                       % Gravitational parameter of the first primary
    mu_r(2) = mu;                                         % Gravitational parameter of the second primary
    R(1:3,:) = [-mu; 0; 0];                               % Location of the first primary in Cartesian coordinates
    R(4:6,:) = [1-mu; 0; 0];                              % Location of the second primary in Cartesian coordinates
    R = repmat(R, 1, size(s,2));                          % Appropriate sizing for the complete state evolutionç
    R = R-repmat(St(1:3,:),2,1);                          % Non-autonomous primaries definition
    
    % Transformation to Cartesian coordinates
    s = cylindrical2cartesian(s, true);
    
    Rr(1:3,:) = s(1:3,:)-R(1:3,:);                        % Relative position to the first primary
    Rr(4:6,:) = s(1:3,:)-R(4:6,:);                        % Relative position to the first primary

    % First primary gravitational acceleration
    q = -dot(s(1:3,:)-2*R(1:3,:),s(1:3,:),1)./sqrt(dot(Rr(1:3,:),Rr(1:3,:),1)).^2;
    f = q.*(3*(1+q)+q.^2)./(1+(1+q).^(3/2));
    gamma(1:3,:) = -mu_r(1)./sqrt(dot(R(1:3,:),R(1:3,:),1)).^3.*((1+f).*s(1:3,:)-f.*R(1:3,:));

    % Second primary gravitational acceleration
    q = -dot(s(1:3,:)-2*R(4:6,:),s(1:3,:),1)./sqrt(dot(Rr(4:6,:),Rr(4:6,:),1)).^2;
    f = q.*(3*(1+q)+q.^2)./(1+(1+q).^(3/2));
    gamma(4:6,:) = -mu_r(2)./sqrt(dot(R(4:6,:),R(4:6,:),1)).^3.*((1+f).*s(1:3,:)-f.*R(4:6,:));

    % Gyroscopic acceleration
    omega = [-2*s(5,:)-s(1,:); 2*s(4,:)-s(2,:); zeros(1,size(s,2))];
    
    % Final vectorfield
    f = gamma(1:3,:)+gamma(4:6,:);                       % Forces field
    a = s(7:9,:)+omega;                                  % Inertial acceleration field

    % Compute the control vector as a dynamics residual
    u = a-f;
end