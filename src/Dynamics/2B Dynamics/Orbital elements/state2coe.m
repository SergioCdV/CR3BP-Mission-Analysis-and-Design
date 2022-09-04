%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/10/21
% File: state2coe.m 
% Issue: 0 
% Validated: 

%% State vector to classical orbital elements %%
% This file contains the function to change from the state vector to the classical orbital elements

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - vector s, the inertial state vector to be converted (position and velocity)
%         - string frame, indicating in which reference frame the state
%           vector is expressed (inertial or perifocal)

% Ouputs: - vector elements, containing the mean classical Euler orbital elements (a, e, RAAN, i, omega, M)

% Units are in S.I

% New version updates:

function [elements] = state2coe(mu, s, frame)
    % Main computation 
    switch (frame)
        case 'Inertial'
            elements = rv2coe(mu,s);
        case 'Perifocal' 
            elements  = prv2coe(mu,s);
        otherwise
            error('No valid reference frame is used');
    end
end

%% Auxiliary functions 
% Function to compute the orbital elements from the inertial state vector 
function [elements] = rv2coe(mu, s)
    % State variables 
    r = s(1:3);                                     % Position vector
    v = s(4:6);                                     % Velocity vector 
    
    % Compute the eccentricy and angular momentum vectors 
    h = cross(r,v);                                 % Angular momentum vector
    e = cross(v,h)/mu-r/norm(r);                    % Eccentricity vector
    K = [0; 0; 1];                                  % Inertial Z axis unit vector
    n = cross(K, h);                                % Node vector
    
    % Compute orbital energy 
    H = norm(v)^2/2-mu/norm(r);
    
    % Determine the type of orbit 
    if (norm(e) ~= 1)
        a = -mu/(2*H);                              % Semimajor axis of the orbit
        p = a*(1-norm(e)^2);                        % Semilatus rectum of the orbit
    else
        p = norm(h)^2/mu;                           % Semilatus rectum of the orbit
        a = Inf;                                    % Semimajor axis of the orbit
    end
    
    % Compute the unit perifocal triad 
    i = e/norm(e); 
    k = h/norm(h); 
    j = cross(k,i);
    
    %Compute the rotation matrix
    Q = [i.'; j.'; k.'];                            % Perifocal rotation matrix
    r0 = Q*r;                                       % Position in the perifocal frame 
    
    %Compute the rest of the orbital elements
    RAAN = atan2(Q(3,1),-Q(3,2));                   % RAAN
    omega = atan2(Q(1,3),Q(2,3));                   % Argument of perigee
    i = acos(Q(3,3));                               % Inclination of the orbit
    
    % Mean anomaly
    theta = atan2(r0(2), r0(1));                    % True anomaly of the orbit

    if (norm(e) >= 1)
        M = norm(e)*sinh(H)-H;                                               % Mean anomaly
    else
        sinE = sqrt(1-norm(e)^2)*sin(theta)/(1+norm(e)*cos(theta));          % Sine of the eccentric anomaly
        cosE = (norm(e)+cos(theta))/(1+norm(e)*cos(theta));                  % Cosine of the eccentric anomaly
        E = atan2(sinE, cosE);                                               % Eccentric anomaly
        M = E-norm(e)*sin(E);                                                % Mean anomaly
    end
        
    % Save the classical orbital elements 
    elements = [a norm(e) RAAN i omega M p];
    elements = rv_singularity(e, n, r, Q, elements);
end

% Function to compute the orbital elements from the perifocal state vector
function [elements] = prv2coe(mu, s)
    % State variables 
    r = s(1:3).';                   % Position vector
    v = s(4:6).';                   % Velocity vector 
    
    % Compute the eccentricy and angular momentum vectors 
    h = cross(r,v);                 % Angular momentum vector
    h = norm(h);                    % Angular momentum norm
    
    % Compute the true anomaly 
    theta = atan2(r(2), r(1)); 
    
    % Compute the eccentricity 
    e = v(2)*(h/mu)-cos(theta); 
    
    % Compute orbital energy 
    H = norm(v)^2/2-mu/norm(r);
    
    % Determine type of orbit 
    if (e ~= 1)
        a = -mu/(2*H);              % Semimajor axis of the orbit
        p = a*(1-e^2);              % Semilatus rectum of the orbit
    else
        p = h^2/mu;                 % Semilatus rectum of the orbit
        a = Inf;                    % Semimajor axis of the orbit
    end
        
    % Compute the rest of the orbital elements
    RAAN = NaN;                     % RAAN
    omega = NaN;                    % Argument of perigee
    i = NaN;                        % Inclination of the orbit
    
    % Mean anomaly
    sinE = sqrt(1-e^2)*sin(theta)/(1+e*cos(theta));     % Sine of the eccentric anomaly
    cosE = (e+cos(theta))/(1+e*cos(theta));             % Cosine of the eccentric anomaly
    E = atan2(sinE, cosE);                              % Eccentric anomaly
    M = E-norm(e)*sin(E);                               % Mean anomaly
    
    % Save the classical orbital elements 
    elements = [a e RAAN i omega M p];
    
    % Singularity warnings 
    tol = 1e-10;               % Circular orbit tolerance    
    if (abs(e) < tol)
        warning('Orbit is circular to numerical precision');   
    end
end

% Function to handle orbital elements singularities when converting from the inertial state vector
function [elements] = rv_singularity(ev, n, r, Q, elements)
    % State variables 
    e = elements(2);           % Orbit eccentricity
    i = elements(4);           % Orbit inclination 
    M = elements(6);           % Mean anomaly
    tol = 1e-10;               % Circular orbit tolerance
    
    % Singularity warnings 
    if (any(isnan(Q)))
        warning('Euler angles are numerically ill-conditioned');
    end
    
    if (abs(e) < tol)
        warning('Orbit is circular to numerical precision');
    end
    
    if (abs(i) < tol)
        warning('Orbit is equatorial to numerical precision');
    end
    
    % Singularity handling 
    if (abs(e) < tol)
        if (abs(i) < tol)
            M = acos(r(1)/norm(r));                 % Circular equatorial orbit
            if (r(2) < 0)
                M = 2*pi-M;
            end
        else
            M = acos(dot(n,r)/(norm(n)*norm(r)));   % Circular inclined orbit
            if (r(3) < 0)
                M = 2*pi-M;
            end
        end
    else
        if (abs(i) < tol)
            M = acos(ev(1)/norm(ev));               % Equatorial orbit
            if (ev(2) < 0)
                M = 2*pi-M; 
            end
        end
    end
    
    % Reconstruction of the orbital elements 
    elements(6) = M;
end