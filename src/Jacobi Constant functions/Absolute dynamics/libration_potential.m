%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 09/02/22
% File: libration_potential.m 
% Issue: 0 
% Validated: 

%% Libration potential %%
% For a given gravitational parameter mu and position vector r, this function computes the 
% Lagrangian function associated to a given libration point. 

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - vector L, the libration point and associated chracteristic
%           distance of the system to be considered.
%         - vector s, a 12x1 array containing the relative synodic position 
%           velocity vectors of the chaser and the target.
%         - scalar order, the order of the approximation in the Legendre
%           series

% Outputs: - scalar L, the Lagrangian function associated with the input state vector. 

% New versions:

function [L] = libration_potential(mu, L, s, order)
    %Constants of the problem 
    mu_r(1) = 1-mu;                     %Gravitational parameters of the most massive primary
    mu_r(2) = mu;                       %Gravitational parameters of the least massive primary
    
    %Location of the unsteady primaries 
    R(1:3,1) = [-mu; 0; 0];             %Location of the first primary
    R(1:3,2) = [1-mu; 0; 0];            %Location of the second primary

    %State variables 
    r = s(1:3);                         %Position vector
    v = s(4:6);                         %Velocity vector

    %Relative state variables 
    U = 0;                              %Relative potential energy
    LID = L(1);                         %Identifier of the libration point
    switch (LID)
        case 1
            gamma = L(end);             %Characteristic distance of the libration point
            r(1) = r(1)-1+mu+gamma;     %Relative x coordinate to the libration point
            r = r/gamma;                %Scaled position vector
            c = legendre_coefficients(mu, LID, gamma, order);

            %Linear potential energy
            for i = 2:order
                ci = c(i);
                P = legendre_polynomials(i+1,r(1)/norm(r));
                U = U + ci*norm(r)^i*P(end);
            end

        case 2
            gamma = L(end);             %Characteristic distance of the libration point
            r(1) = r(1)-1+mu-gamma;     %Relative x coordinate to the libration point
            r = r/gamma;                %Scaled position vector
            c = legendre_coefficients(mu, LID, gamma, order);

            %Potential energy
            for i = 1:order
                ci = c(i);
                P = legendre_polynomials(i+1,r(1)/norm(r));
                U = U + ci*norm(r)^i*P(end);
            end

        case 3
            gamma = L(end);             %Characteristic distance of the libration point
            r(1) = r(1)+mu+gamma;       %Relative x coordinate to the libration point
            r = r/gamma;                %Scaled position vector

            %Potential energy
            for i = 2:order
                c = legendre_coefficients(mu, LID, gamma, i);
                c = c(end);
                P = legendre_polynomials(i+1,r(1)/norm(r));
                U = U + c*norm(r)^i*P(end);
            end

        case 4
            l = [mu+1/2; sqrt(3)/2; 0];     %Libration point position vector
            r = r - l;                      %Relative x coordinate to the libration point
            R = R - repmat(l,1,2);          %Relative primaries position vectors

           %Potential energy
            for i = 1:length(mu_r)
                for j = 1:order+1
                    P = legendre_polynomials(j,dot(R(:,i),r)/(norm(R(:,i)*norm(r))));
                    U = U + (norm(r)/norm(R(:,i)))^(j-1)*P(end)/norm(R(:,i));
                end
            end

        case 5
            l = [mu+1/2; -sqrt(3)/2; 0];    %Libration point position vector
            r = r - l;                      %Relative x coordinate to the libration point
            R = R - repmat(l,1,2);          %Relative primaries position vectors

            %Potential energy
            for i = 1:length(mu_r)
                for j = 1:order+1
                    P = legendre_polynomials(j,dot(R(:,i),r)/(norm(R(:,i)*norm(r))));
                    U = U + (norm(r)/norm(R(:,i)))^(j-1)*P(end)/norm(R(:,i));
                end
            end

        otherwise
            error('Libration point calculation has not been implemented yet');
    end

    %Total kinetic energy
    v = v + [-r(2); r(1); 0];           %Total velocity vector    
    T = (1/2)*dot(v,v);

    %Relative potential function 
    L = T+U;
end