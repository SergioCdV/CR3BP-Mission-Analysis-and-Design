%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 08/05/21
% File: LQR_control.m 
% Issue: 0 
% Validated: 08/05/21

%% State Dependt Ricatti Equation Control %%
% This script contains the function to compute the control law by means of an SDRE controller.

% Inputs: - string model, selecting the linear model to compute the linear state
%           transition matrix of the system
%         - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - array Sg, the guidance law to follow
%         - array Sn, the system state
%         - array St, the target orbit reference state 
%         - scalar Ln, the libration point number. It may be left 0 if the
%           target orbit is not librating around any Lagrange point
%         - scalar gamma, the relative distance of the Ln point to the
%           nearest primary. Again, it may be left as 0 if needed
%         - matrices Q and M, penalizing on the state error and the control
%           effort

% Output: - vector u, the computed control law

% New versions: 

function [u] = SDRE_control(model, mu, Sg, Sn, St, Ln, gamma, Q, M)
    % Approximation 
    n = 6;                              % Dimension of the approximation
    order = 2;                          % Order of the approximation 

    % Model coefficients and constants
    mup(1) = 1-mu;                      % Reduced gravitational parameter of the first primary 
    mup(2) = mu;                        % Reduced gravitational parameter of the second primary 
    R(:,1) = [-mu; 0; 0];               % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];              % Synodic position of the second primary
    
    % Linear model matrices
    B = [zeros(n/2); eye(n/2); zeros(n/2)];         % Linear model input matrix 
    Omega = [0 2 0; -2 0 0; 0 0 0];                 % Coriolis dyadic
    
    % Sanity check on the libration point 
    if (Ln == 0)
        if (gamma == 0)
            model = 'RLM'; 
        else
            error('No valid linear model was selected');
        end
    end 
    
    % Preallocation 
    u = zeros(3,size(Sn,1));                        % Control law

    % Compute the trajectory
    for i = 1:size(Sn,1)
        % State coefficients 
        r_t = St(i,1:3).';                          % Synodic position of the target

        % Select linear model 
        switch (model)
            case 'SLLM'
                cn = legendre_coefficients(mu, Ln, gamma, order);     % Compute the relative Legendre coefficient c2 
                c2 = cn(2); 
                Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              % Linear model state matrix

                % Linear state model
                A = [zeros(3) eye(3) zeros(3); Sigma Omega zeros(3); eye(3) zeros(3) zeros(3)]; 
                
            case 'ULLM' 
                cn = relegendre_coefficients(mu, r_t.', order);       % Compute the relative Legendre coefficient c2 
                c2 = cn(2); 
                Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              % Linear model state matrix

                % Linear state model
                A = [zeros(3) eye(3) zeros(3); Sigma Omega zeros(3); eye(3) zeros(3) zeros(3)]; 
                
            case 'RLM' 
                % Relative position between the primaries and the target 
                Ur(:,1) = r_t-R(:,1);               % Position of the target with respect to the first primary
                ur(:,1) = Ur(:,1)/norm(Ur(:,1));    % Unit vector of the relative position of the target with respect to the primary
                Ur(:,2) = r_t-R(:,2);               % Position of the target with respect to the first primary
                ur(:,2) = Ur(:,2)/norm(Ur(:,2));    % Unit vector of the relative position of the target with respect to the primary

                % Evaluate the linear model 
                Sigma = -((mup(1)/norm(Ur(:,1))^3)+(mup(2)/norm(Ur(:,2)))^3)*eye(3) ...
                        +3*((mup(1)/norm(Ur(:,1))^3)*(ur(:,1)*ur(:,1).')+(mup(2)/norm(Ur(:,2))^3)*(ur(:,2)*ur(:,2).'));

                % Linear state model
                A = [zeros(3) eye(3) zeros(3); Sigma Omega zeros(3); eye(3) zeros(3) zeros(3)]; 

            case 'Numerical'
                % Linear state model
                A = [rel_jacobian(mu, [St Sn].') zeros(6,3); eye(3) zeros(3) zeros(3)];  

            otherwise 
                error('No valid linear model was selected'); 
        end

        % Compute the feedback control law
        [K,~,~] = lqr(A,B,Q,M);
        u(:,i) = -K*(Sn(i,:)-Sg(i,:)).';  
    end
end