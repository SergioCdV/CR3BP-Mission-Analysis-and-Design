%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/04/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the 9xm state vector 
%         - scalar tf, the final time of flight
%         - string method, the collocation distribution to be used

% Outputs: - vector u, the nondimensional 3xm control vector

function [u] = acceleration_control(mu, St, C, tf, method)
    % Constants of the problem 
    mu_r(1) = 1-mu;                 % Gravitational parameter of the first primary 
    mu_r(2) = mu;                   % Gravitational parameter of the second primary 
    R(1:3,:) = [-mu; 0; 0];         % Location of the first primary in Cartesian coordinates
    R(4:6,:) = [1-mu; 0; 0];        % Location of the second primary in Cartesian coordinates
    R = repmat(R, 1, size(C,2));    % Appropriate sizing for the complete state evolution
    R = R-[St(1:3,:); St(1:3,:)];   % Non-autonomous primaries definition

    % Transformation to Cartesian coordinates
    s = cylindrical2cartesian(C, true);    

    % Compute the radius vector
    r = sqrt(s(1,:).^2+s(2,:).^2+s(3,:).^2);

    % Compute the control vector as a residual of the dynamics
    switch (method)
        case 'Sundman'
            % Normalizing factor
            c = tf;

            % Derivative of the radius with the arclength 
            dr = dot(C(1,:).*C(4,:)+C(3,:).*C(6,:))./r;

            % Compute the control vector as a dynamics residual
            u = [C(7,:)-dr.*C(4,:)./r+c.^2.*mu.*C(1,:)./r-C(1,:).*C(5,:).^2; ...
                 C(1,:).*C(8,:)+2*C(4,:).*C(5,:); ... 
                 C(9,:)-dr.*C(6,:)./r+c.^2.*mu.*C(3,:)./r];

            u(2,:) = u(2,:)./r.^2; 

        otherwise
            % Normalizing factor
            c = tf;

            % Compute the control vector as a dynamics residual of the Encke vector field
            Rr(1:3,:) = s(1:3,:)-R(1:3,:);                                                  % Relative position to the first primary
            Rr(4:6,:) = s(1:3,:)-R(4:6,:);                                                  % Relative position to the first primary

            % First primary gravitational acceleration
            q = -dot(2*Rr(1:3,:),s(1:3,:))./sqrt(dot(Rr(1:3,:),Rr(1:3,:),2)).^2;
            f = q.*(3*(1+q)+q.^2)./(1+(1+q).^(3/2));
            gamma(1:3,:) = -mu_r(1)/sqrt(dot(R(1:3,:),R(1:3,:),2)).^3*((1+f)*s(1:3,:)-f*R(1:3,:)); 

            % Second primary gravitational acceleration
            q = -dot(2*Rr(4:6,:),s(1:3,:))./sqrt(dot(Rr(4:6,:),Rr(4:6,:),2)).^2;
            f = q.*(3*(1+q)+q.^2)./(1+(1+q).^(3/2));
            gamma(4:6,:) = -mu_r(2)/sqrt(dot(R(4:6,:),R(4:6,:),2)).^3*((1+f)*s(1:3,:)-f*R(4:6,:));  

            % Gyroscopic acceleration
            omega = [-2*s(5,:)-s(1,:); 2*s(4,:)-s(2,:); zeros(1,size(C,2))]; 

            % Final vectorfield
            u = s(7:9,:)+omega+c^2*(gamma(1:3,:) + gamma(4:6,:));                           
    end
end
