%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 09/02/22
% File: libration_potential.m 
% Issue: 0 
% Validated: 

%% Libration potential %%
% For a given gravitational parameter mu and position vector r, this function computes the 
% Lagrangian function associated to a given libration point, together with
% the associated high order potential.

% Inputs: - scalar mu, the reduced gravitational parameter of the system. 
%         - vector L, the libration point and associated chracteristic
%           distance of the system to be considered.
%         - vector s, a 12x1 array containing the relative synodic position 
%           velocity vectors of the chaser and the target.
%         - scalar order, the order of the approximation in the Legendre
%           series

% Outputs: - scalar L, the Lagrangian function associated with the input state vector.
%          - scalar U, the potential function associated to the libration
%            point

% New versions:

function [L, U, Uh] = libration_potential(mu, L, s, order)
    %Constants of the problem 
    mu_r(1) = 1-mu;                     %Gravitational parameters of the most massive primary
    mu_r(2) = mu;                       %Gravitational parameters of the least massive primary
    
    %Location of the unsteady primaries 
    R(1:3,1) = [-mu; 0; 0];             %Location of the first primary
    R(1:3,2) = [1-mu; 0; 0];            %Location of the second primary

    %State variables 
    r = s(1:3);                         %Position vector
    v = s(4:6);                         %Velocity vector

    %High order potential energy 
    Uh = [0; 0];                        %High order potential initialization
    rmag = norm(r);                     %Norm of the position vector
    Rmag(1) = norm(R(:,1));             %Norm of the first primary position vector
    Rmag(2) = norm(R(:,2));             %Norm of the second primary position vector

    for i = 1:length(mu_r)
        if (norm(r)/Rmag(i) < 1)
            for j = 2:order
                P = legendre_polynomials(j+1, dot(r,R(:,i))/(rmag*Rmag(i)));
                Uh(i) = Uh(i) + (norm(r)/Rmag(i))^(j)*P(end);
            end
            Uh(i) = mu_r(i)*Uh(i)/Rmag(i);
        else
            Uh(i) = mu_r(i)/norm(r-R(:,i));
        end
    end
    Uh = -sum(Uh);

    %Relative state variables 
    LID = L(1);                                         %Identifier of the libration point
    gamma = L(end);                                     %Characteristic distance of the libration point
    r = synodic2lagrange(mu, gamma, LID, r, true);      %Normalized relative position vector to the libration point

    switch (LID)
        case 1
            %Potential energy
            c = legendre_coefficients(mu, LID, gamma, 2);               %Characteristic frequencies
            c = c(end);                                                 %Characteristic c2 frequency
            U = -c*r(1)^2+(c/2)*(r(2)^2+r(3)^2);                        %Potential energy

        case 2
            %Potential energy
            c = legendre_coefficients(mu, LID, gamma, 2);               %Characteristic frequencies
            c = c(end);                                                 %Characteristic c2 frequency
            U = -c*r(1)^2+(c/2)*(r(2)^2+r(3)^2);                        %Potential energy

        case 3
            %Potential energy
            c = legendre_coefficients(mu, LID, gamma, 2);               %Characteristic frequencies
            c = c(end);                                                 %Characteristic c2 frequency
            U = -c*r(1)^2+(c/2)*(r(2)^2+r(3)^2);                        %Potential energy

        case 4
            %Potential energy
            a = 3*sqrt(3)/4*(1-2*mu);                                   %Invariant manifolds constant
            U = (1/8)*r(1)^2-(5/8)*r(2)^2+(1/2)*r(3)^2-a*r(1)*r(2);     %Linear potential energy

        case 5
            %Potential energy
            a = -3*sqrt(3)/4*(1-2*mu);                                  %Invariant manifolds constant
            U = (1/8)*r(1)^2-(5/8)*r(2)^2+(1/2)*r(3)^2-a*r(1)*r(2);     %Linear potential energy

        otherwise
            error('No valid libration point was selected');
    end

    %Total kinetic energy
    v = v + [-r(2); r(1); 0];           %Total velocity vector    
    T = (1/2)*dot(v,v);

    %Relative potential function 
    L = T-U;
end