%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/07/21
% File: controlabillity.m 
% Issue: 0 
% Validated: 24/07/21

%% Controlabillity %%
% This script contains the function to study the controlabillity of a given linear model

% Inputs: - string model, the linear dynamics to be analysed
%         - scalar mu, the reduced gravitational parameter of the system
%         - array Sn, containing the target spacecraft trajectory along
%           which to evaluate the model
%         - scalar index, the number of phase space points to be analysed
%         - scalar Ln, indicating the associated target's libration point
%         - scalar gamma, the distance to the nearest primary from Ln

% Output: - vector contrabillity, a boolean array to indicate whether the
%           model is controlable or not phase space point-wise

% New versions: 

function [controlable] = controlabillity(model, mu, Sn, index, Ln, gamma)
    % Approximation 
    n = 6;                              % Dimension of the state vector
    order = 2;                          % Order of the approximation 

    % Preallocation 
    controlable = zeros(index,1);       % Controllability gramian

    % Model coefficients 
    mup(1) = 1-mu;                      % Reduced gravitational parameter of the first primary 
    mup(2) = mu;                        % Reduced gravitational parameter of the second primary 
    R(:,1) = [-mu; 0; 0];               % Synodic position of the first primary
    R(:,2) = [1-mu; 0; 0];              % Synodic position of the second primary

    % Linear model matrices
    B = [zeros(n/2); zeros(n/2); eye(n/2)];         % Linear model input matrix 
    Omega = [0 2 0; -2 0 0; 0 0 0];                 % Coriolis dyadic

    for i = 1:index
        % State coefficients 
        r_t = Sn(i,1:3).';                          % Synodic position of the target
        
        % Select linear model 
        switch (model)
            case 'SLLM'
                cn = legendre_coefficients(mu, Ln, gamma, order);     % Compute the absolute Legendre coefficient c2 
                c2 = cn(2); 
                Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              % Linear model state matrix
            case 'ULLM' 
                cn = relegendre_coefficients(mu, r_t.', order);       % Compute the relative Legendre coefficient c2 
                c2 = cn(2); 
                Sigma = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];              % Linear model state matrix
            case 'RLM' 
                % Relative position between the primaries and the target 
                Ur1 = r_t-R(:,1);               % Position of the target with respect to the first primary
                ur1 = Ur1/norm(Ur1);            % Unit vector of the relative position of the target with respect to the primary
                Ur2 = r_t-R(:,2);               % Position of the target with respect to the first primary
                ur2 = Ur2/norm(Ur2);            % Unit vector of the relative position of the target with respect to the primary

                %Evaluate the linear model 
                Sigma = -((mup(1)/norm(Ur1)^3)+(mup(2)/norm(Ur2))^3)*eye(3)+3*((mup(1)/norm(Ur1)^3)*(ur1*ur1.')+(mup(2)/norm(Ur2)^3)*(ur2*ur2.'));

            otherwise 
                error('No valid linear model was selected'); 
        end

        % Linear state model
        A = [zeros(3) eye(3) zeros(3); zeros(3) zeros(3) eye(3); zeros(3) Sigma Omega];     

        % Controlability matrix 
        C = ctrb(A,B); 
        controlable(i) = (rank(C) == size(A,1));
    end
end