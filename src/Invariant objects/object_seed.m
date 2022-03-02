%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/20
% File: object_seed.m 
% Issue: 0 
% Validated: 

%% Linear seed %%
% This script provides a function to generate a linear seed for various 
% objects or dynamical solutions.

% Inputs: - scalar mu, the reduced gravitational parameter of the system.
%         - vector parameters, containing the solution characteristic
%           variables, as specified at each seed case.
%         - string object, to select the seed case to apply. Currently only
%           halo and quasi-periodic torus seeds are available.

% Output: - vector field seed, containing the required initial solution
%           seed. 

% Methods: - Lyapunov orbits seeds employ a first order approximation from
%           (solutions of the linearized equations of motion around the Lagrange
%            point Li) from Howell et al.
%          - Halo orbits seeds employ a semianalytical 3th order
%            Richardson approximation. 

% New versions: include L3 solutions.

function [seed, T] = object_seed(mu, parameters, object)
    %Main interface 
    switch (object) 
        case 'Lyapunov' 
            [seed, T] = lyapunov_seed(mu, parameters);         %Generate a Lyapunov linear seed
        case 'Halo' 
            [seed, T] = halo_seed(mu, parameters);             %Generate a Halo 3D third order seed            
        otherwise 
            error('No valid options was selected');            %Display selection error
    end
end

%% Auxiliary functions 
%Lyapunov linear seed orbit 
function [seed, T] = lyapunov_seed(mu, parameters)
    %Constants 
    rho = 10000;            %Number of points per period
    
    %Parameters of the Lyapunov orbit
    Ax = parameters(1);     %In-plane trajectory
    Az = parameters(2);     %Out-of-plane trajectory     
    phi = parameters(3);    %In-plane phase
    psi = parameters(4);    %Out-of-plane phase
    L = parameters(5);      %Lagrange point identifier
    gamma = parameters(6);  %Lagrange point coordinate
    n = parameters(7);      %Number of periods to generate
        
    %Orbit parameters (frequencies)
    cn = legendre_coefficients(mu, L, gamma, 2);                %Legendre coefficient c_2 (equivalent to mu)
    c2 = cn(2);                                                 %Legendre coefficient c_2 (equivalent to mu)
    wp  = sqrt((1/2)*(2-c2+sqrt(9*c2^2-8*c2)));                 %In-plane frequency
    wv  = sqrt(c2);                                             %Out of plane frequency
    kap = (wp^2+1+2*c2)/(2*wp);                                 %Contraint on the planar amplitude
    
    %Temporal parametrization
    T = (2*pi)/wp;                     %Period of the orbit 
    tspan = linspace(0, n*T, rho);     %Integration vector
    
    %Seed trajectory
    x = -Ax*cos(wp*tspan+phi);         %X relative coordinate
    y = kap*Ax*sin(wp*tspan+phi);      %Y relative coordinate
    z = Az*sin(wv*tspan+psi);          %Z relative coordinate
    vx = wp*Ax*sin(wp*tspan+phi);      %Vx relative velocity
    vy = kap*wp*Ax*cos(wp*tspan+phi);  %Vy relative velocity
    vz = wv*Az*cos(wv*tspan+psi);      %Vz relative velocity  
    
    %Relative to synodic reference frame transformation
    if (L == 1) 
        k = -1;
    elseif (L == 2)
        k = 1;
    else
        error('No valid Lagrange point was selected');
    end
    
    x = gamma*x+(1-mu+k*gamma);             %Synodic X coordinate
    y = gamma*y;                            %Synodic Y coordinate
    z = gamma*z;                            %Synodic Z coordinate
    vx = gamma*vx;                          %Synodic x velocity
    vy = gamma*vy;                          %Synodic y velocity
    vz = gamma*vz;                          %Synodic z velocity 
    
    %Output seed
    seed = [x.' y.' z.' vx.' vy.' vz.'];    %Seed trajectory
end

%Halo third order seed orbit
function [seed, T] = halo_seed(mu, parameters)
    %Constants 
    dtheta = 1e-3;              %Step in the angular variable
    
    %Parameters of the orbit 
    n = parameters(1);          %+1 for northern halo, -1 for southern halo
    Az = parameters(2);         %Out-of-plane amplitude
    L = parameters(3);          %Lagrange point identifier
    gamma = parameters(4);      %Lagrange point coordinate
    m = parameters(5);          %Number of periods
    
    %Independent variable     
    tau = 0:dtheta:2*pi*m;
    
    %Determine some boolean parameters for the halo determination, concerning the nondimensional reference frame used in the
    %approximation
    if (L == 1) 
        won = 1;            %Associated sign
        primary = 1-mu;     %Reference primary position
    elseif (L == 2) 
        won = -1;           %Associated sign
        primary = 1-mu;     %Reference primary position
    elseif (L == 3) 
        won = 1;            %Associated sign
        primary = -mu;      %Reference primary position
    else
        error('No valid Lagrange point was selected'); 
    end
    
    %Legendre polynomial coefficients c_n for the Richardson approximation
    order = 4;                                          %Order of the approximation
    c = legendre_coefficients(mu, L, gamma, order);     %Legendre coefficients
    
    %Determine the orbit spatial eigenvalue 
    polylambda = [1 0 (c(2)-2) 0 -(c(2)-1)*(1+2*c(2))];
    lambda = sort(roots(polylambda));
    if (L == 3) 
        lambda = abs(lambda(3));
    else
        lambda = abs(lambda(1));
    end

    %Richardson 3th order approximation coefficients
    k = 2*lambda/(lambda^2+1-c(2));
    del = lambda^2-c(2);
    
    d1 = ((3*lambda^2)/k)*(k*(6*lambda^2-1)-2*lambda);
    d2 = ((8*lambda^2)/k)*(k*(11*lambda^2-1)-2*lambda);
    
    a21 = (3*c(3)*(k^2-2))/(4*(1+2*c(2)));
    a22 = 3*c(3)/(4*(1+2*c(2)));
    a23 = -(3*c(3)*lambda/(4*k*d1))*(3*(k^3)*lambda-6*k*(k-lambda)+4);
    a24 = -(3*c(3)*lambda/(4*k*d1))*(2+3*k*lambda);
    
    b21 = -(3*c(3)*lambda/(2*d1))*(3*k*lambda-4);
    b22 = 3*c(3)*lambda/d1;
    d21 = -c(3)/(2*lambda^2);
    
    a31 = -(9*lambda/(4*d2))*(4*c(3)*(k*a23-b21)+k*c(4)*(4+k^2))+((9*lambda^2+1-c(2))/(2*d2))*(3*c(3)*(2*a23-k*b21)+c(4)*(2+3*k^2));
    a32 = -(1/d2)*((9*lambda/4)*(4*c(3)*(k*a24-b22)+ k*c(4))+1.5*(9*lambda^2+1-c(2))*(c(3)*(k*b22+d21-2*a24)-c(4)));
    
    b31 = (0.375/d2)*(8*lambda*(3*c(3)*(k*b21-2*a23)-c(4)*(2+3*k^2))+(9*lambda^2+1+2*c(2))*(4*c(3)*(k*a23-b21)+k*c(4)*(4+ k^2)));
    b32 = (1/d2)*(9*lambda*(c(3)*(k*b22+d21-2*a24)-c(4))+0.375*(9*lambda^2+1+2*c(2))*(4*c(3)*(k*a24-b22)+k*c(4)));
    
    d31 = (3/(64*lambda^2))*(4*c(3)*a24+c(4));
    d32 = (3/(64*lambda^2))*(4*c(3)*(a23-d21)+c(4)*(4+k^2));
    
    s1 = (1/(2*lambda*(lambda*(1+k^2)-2*k)))*(1.5*c(3)*(2*a21*(k^2-2)-a23*(k^2+2)-2*k*b21)-0.375*c(4)*(3*k^4-8*k^2+8));
    s2 = (1/(2*lambda*(lambda*(1+k^2)-2*k)))*(1.5*c(3)*(2*a22*(k^2-2)+a24*(k^2+2)+2*k*b22+5*d21)+0.375*c(4)*(12-k^2));
    
    a1 = -1.5*c(3)*(2*a21+a23+5*d21)-0.375*c(4)*(12-k^2);
    a2 = 1.5*c(3)*(a24-2*a22)+1.125*c(4);
    l1 = a1 + 2*(lambda^2)*s1;
    l2 = a2 + 2*(lambda^2)*s2;
    deltan = won*n;

    %In-plane amplitude (related to Az by a non-linear analytical constraint)
    Ax = sqrt((-del-l2*Az^2)/l1);

    %Phase space vector
    x = a21*Ax^2+a22*Az^2-Ax*cos(tau)+(a23*Ax^2-a24*Az^2)*cos(2*tau)+(a31*Ax^3-a32*Ax*Az^2)*cos(3*tau);
    y = k*Ax*sin(tau)+(b21*Ax^2-b22*Az^2)*sin(2*tau)+(b31*Ax^3-b32*Ax*Az^2)*sin(3*tau);
    z = deltan*Az*cos(tau)+deltan*d21*Ax*Az*(cos(2*tau)-3)+deltan*(d32*Az*Ax^2-d31*Az^3)*cos(3*tau);
    dx = lambda*Ax*sin(tau)-2*lambda*(a23*Ax^2-a24*Az^2)*sin(2*tau)-3*lambda*(a31*Ax^3-a32*Ax*Az^2)*sin(3*tau);
    dy = lambda*(k*Ax*cos(tau)+2*(b21*Ax^2-b22*Az^2)*cos(2*tau)+3*(b31*Ax^3-b32*Ax*Az^2)*cos(3*tau));
    dz = -lambda*deltan*Az*sin(tau)-2*lambda*deltan*d21*Ax*Az*sin(2*tau)-3*lambda*deltan*(d32*Az*Ax^2-d31*Az^3)*sin(3*tau);

    %Position vector
    x = (gamma*(x-won)+primary);   %Rescaled synodic x coordinate
    y = gamma*y;                   %Rescaled synodic y coordinate
    z = gamma*z;                   %Rescaled synodic z coordinate
    r0 = [x.' -y.' z.'];           %Position vector

    %Velocity vector
    v0 = gamma*[dx.' dy.' dz.'];                            
    
    %Output
    seed = [r0 v0];
    T = 2*pi/lambda;
end