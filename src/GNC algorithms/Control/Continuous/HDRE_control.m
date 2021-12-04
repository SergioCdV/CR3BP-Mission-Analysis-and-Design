%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 04/12/21
% File: HDRE_control.m 
% Issue: 0 
% Validated: 04/12/21

%% Hinfty SDRE Control %%
% This script contains the function to compute the control law by means of an SDRE-H controller.

% Inputs: - string model, selecting the linear model to compute the linear state
%           transition matrix of the system
%         - scalar mu, the reduced gravitational parameter of the CR3BP
%           system
%         - array St, the target orbit reference state 
%         - array Sg, the guidance law to follow
%         - array Sn, the system state
%         - scalar Ln, the libration point number. It may be left 0 if the
%           target orbit is not librating around any Lagrange point
%         - scalar gamma, the relative distance of the Ln point to the
%           nearest primary. Again, it may be left as 0 if needed
%         - boolean SDRE, for time-varying dynamics
%         - vector W, variance of the noise/disturbance vector

% Output: - vector u, the computed control law

% New versions: 

function [u] = HDRE_control(model, mu, Sg, Sn, St, Ln, gamma, SDRE, W)
   %Approximation 
    n = 6;                              %Dimension of the approximation
    
    %Sanity check on the libration point 
    if ((Ln == 0) && (model ~= 'RLM'))
        if (gamma == 0)
            model = 'RLM'; 
        else
            error('No valid linear model was selected');
        end
    end 
    
    %Preallocation 
    u = zeros(3,size(Sn,1));                        %Control law

    %Define the performance system  
    C = eye(9);                                     %Observer matrix
    D = [zeros(9,6) [zeros(3); eye(3); zeros(3)]];  %Observer control matrix

    %Control input matrix 
    B = [zeros(n/2); eye(n/2); zeros(n/2)];         %Linear model input matrix 
    Bn = [eye(6)*W; zeros(3,6)];                    %Noise input matrix           
    Bt = [Bn 1e6*B];                                %Total input matrix
    ncont = 3;                                      %Number of control inputs

    %Performance 
    pole = 1; 

    %Time-invariant plant
    if (~SDRE)
        %Plant
        A = dynamics(mu, Ln, gamma, St, zeros(6,1), model);
        P = ss(A, Bt, C, D);

        %H-infty controller
        [K, ~, ~] = hinffi(P, ncont, pole);              
    end

    %Compute the trajectory
    for i = 1:size(Sn,1)
        %Time-varying dynamics matrix
        if (SDRE)
            %Plant
            A = dynamics(mu, Ln, gamma, St(i,:), Sn(i,:), model);
            P = ss(A, Bt, C, D);

            %H-infty controller
            [K, ~, ~] = hinffi(P, ncont, pole);              
        end

        %Compute some white noise 
        w = -2*ones(6,1)+rand(6,1);

        %Compute the feedback control law
        u(:,i) = K*[(Sn(i,:)-Sg(i,:)).'; w];      %Compute the control law
    end
end

%% Auxiliary functions 
%Linear system
function [A] = dynamics(mu, Ln, gamma, St, Sn, model)
    %Order of the approximation 
    order = 2;                          

    %Model coefficients 
    mup(1) = 1-mu;                      %Reduced gravitational parameter of the first primary 
    mup(2) = mu;                        %Reduced gravitational parameter of the second primary 
    R(:,1) = [-mu; 0; 0];               %Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];              %Synodic position of the second primary    

    %Coriolis dyadic
    Omega = [0 2 0; -2 0 0; 0 0 0];  

    %State coefficients 
    r_t = St(1:3).';                          %Synodic position of the target

    %Select linear model 
    switch (model)
        case 'SLLM'
            cn = legendre_coefficients(mu, Ln, gamma, order);     %Compute the relative Legendre coefficient c2 
            c2 = cn(2); 
            Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              %Linear model state matrix

            %Linear state model
            A = [zeros(3) eye(3) zeros(3); Sigma Omega zeros(3); eye(3) zeros(3) zeros(3)]; 
            
        case 'ULLM' 
            cn = relegendre_coefficients(mu, r_t.', order);       %Compute the relative Legendre coefficient c2 
            c2 = cn(2); 
            Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              %Linear model state matrix

            %Linear state model
            A = [zeros(3) eye(3) zeros(3); Sigma Omega zeros(3); eye(3) zeros(3) zeros(3)]; 
            
        case 'RLM' 
            %Relative position between the primaries and the target 
            Ur(:,1) = r_t-R(:,1);               %Position of the target with respect to the first primary
            ur(:,1) = Ur(:,1)/norm(Ur(:,1));    %Unit vector of the relative position of the target with respect to the primary
            Ur(:,2) = r_t-R(:,2);               %Position of the target with respect to the first primary
            ur(:,2) = Ur(:,2)/norm(Ur(:,2));    %Unit vector of the relative position of the target with respect to the primary

            %Evaluate the linear model 
            Sigma = -((mup(1)/norm(Ur(:,1))^3)+(mup(2)/norm(Ur(:,2)))^3)*eye(3) ...
                    +3*((mup(1)/norm(Ur(:,1))^3)*(ur(:,1)*ur(:,1).')+(mup(2)/norm(Ur(:,2))^3)*(ur(:,2)*ur(:,2).'));

            %Linear state model
            A = [zeros(3) eye(3) zeros(3); Sigma Omega zeros(3); eye(3) zeros(3) zeros(3)]; 

        case 'Numerical'
            %Linear state model
            A = [rel_jacobian(mu, [St Sn].') zeros(6,3); eye(3) zeros(3) zeros(3)]; 

        otherwise 
            error('No valid linear model was selected'); 
    end
end