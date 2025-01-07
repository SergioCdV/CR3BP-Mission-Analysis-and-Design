%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 31/12/24
% File: OrbitSeed.m 
% Issue: 0 

%% Halo Orbit seed %% 
% This functions allows to generate a 3rd order Halo orbit seed

% Inputs: - obj Orbit, the Lissajous orbit of interest
%         - vector Amp, the amplitudes of the Lissajous orbit
%         - vector theta[2 x N], the set of phases of the orbit 
%         - double order, order of the seed approximation (1 or 3)
%         - vector freq, the set of frequencies of the orbit 
%         - double kap, the amplitude constraint in the xy plane

% Output: - array seed [6 x N], containing the required initial solution seed

function [seed] = OrbitSeed(obj, Amp, theta, order, freq, kap) 
    % Sanity checks 
    if ( ~exist("freq", "var") )
        freq = obj.OrbitFrequencies;
    end

    if ( ~exist("kap", "var") )
        kap = obj.kap;
    end

    if ( ~exist("theta", "var") )
        theta = zeros(2,1);
    end

    if ( ~exist("order", "var") )
        order = 3;
    end

    if ( length(Amp) < 2 )
        O = zeros( 2-size(Amp,1), size(Amp,2) );
        Amp = [Amp; O];
    end

    if ( order == 1 )
        % Lissajous seed
        seed = src.OrbitFamilies.LissajousOrbit( obj.System, obj.LibrationPoint ).OrbitSeed( Amp, theta, freq, kap );
    
    else
        if ( order ~= 3 )
            warning('The input order of the halo orbit seed is not supported. Generating a 3rd order seed...')
        end

        % 3rd order seed 
        [seed, lambda] = richardson_seed(obj.System.mu, obj.LibrationPoint, obj.System.LP.gamma( obj.LibrationPoint ), obj.Branch, Amp, theta);
    end
end

%% Auxiliary function 
function [seed, lambda] = richardson_seed(mu, L, gamma, branch, Amp, tau)
    % Parameters of the orbit 
    n = branch(1);              % +1 for northern halo, -1 for southern halo
    Az = Amp(1);                % Out-of-plane amplitude

    % Dimensionalising
    Az = Az / gamma; 
    
    % Determine some boolean parameters for the halo determination, concerning the nondimensional reference frame used
    switch (L)
        case 1
            won = 1;            % Associated sign

        case 2 
            won = -1;           % Associated sign

        case 3
            won = 1;            % Associated sign

        otherwise 
            error('No valid Lagrange point was selected'); 
    end
    
    % Legendre polynomial coefficients c_n for the Richardson approximation
    order = 4;                                                                  % Order of the approximation
    cn = src.Systems.CR3BPSystem.LegendreCoefficients(mu, L, gamma, order);     % Legendre coefficients
    
    % Determine the orbit spatial eigenvalue
    polylambda = [1 0 (cn(3) - 2) 0 -(cn(3) - 1) * (1 + 2 * cn(3))];
    lambda = roots( polylambda );

    if ( L == 3 )
        lambda = abs( lambda(3) ) ;
    else        
        lambda = abs( lambda(1) ) ;
    end

    % Richardson 3rd order approximation coefficients
    k = 2 * lambda / (lambda^2 + 1 - cn(3));
    del = lambda^2 - cn(3);
    
    d1 = (3 * lambda^2 / k) * (k * (6 * lambda^2 - 1) - 2 * lambda);
    d2 = (8 * lambda^2 / k) * (k * (11 * lambda^2 - 1) - 2 * lambda);
    
    a21 = 3 * cn(4) * (k^2 - 2) / ( 4 * (1 + 2 * cn(3)) );
    a22 = 3 * cn(4) / (4 * (1 + 2 * cn(3)));
    a23 = -(3 * cn(4) * lambda / (4 * k * d1)) * (3 * k^3 * lambda - 6 * k * (k - lambda) + 4);
    a24 = -(3 * cn(4) * lambda / (4 * k * d1)) * (2 + 3 * k * lambda);
    
    b21 = -3 * cn(4) * lambda / (2 * d1) * (3 * k * lambda - 4);
    b22 = 3 * cn(4) * lambda / d1;
    d21 = -cn(4) / (2 * lambda^2);
    
    a31 = -9 * lambda / (4 * d2) * (4 * cn(4) * (k * a23 - b21) + k * cn(5) * (4 + k^2)) + ((9 * lambda^2 + 1 - cn(3)) / (2 * d2)) * (3 * cn(4) * (2 * a23 - k * b21) + cn(5) * (2 + 3 * k^2));
    a32 = -9 * lambda / (4 * d2) * (4 * cn(4) * (k * a24 - b22) + k * cn(5)) - 1.5 * (9 * lambda^2 + 1 - cn(3)) * (cn(4) * (k * b22 + d21 - 2 * a24) - cn(5));
    
    b31 = (0.375 / d2) * (8 * lambda * (3 * cn(4) * (k * b21 - 2 * a23) - cn(5) * (2 + 3 * k^2)) + (9 * lambda^2 + 1 + 2 * cn(3)) * (4 * cn(4) * (k * a23 - b21) + k * cn(5) * (4 + k^2)));
    b32 = 9 * lambda / d2 * (cn(4) * (k * b22 + d21 - 2 * a24) - cn(5)) + 0.375 / d2 * (9 * lambda^2 + 1 + 2 * cn(3)) * (4 * cn(4) * (k * a24 - b22) + k * cn(5));
    
    d31 = (3 / (64 * lambda^2)) * (4 * cn(4) * a24 + cn(5));
    d32 = (3 / (64 * lambda^2)) * (4 * cn(4) * (a23 - d21) + cn(5) * (4 + k^2));
    
    s1 = (1.5 * cn(4) * (2 * a21 * (k^2 - 2) - a23 * (k^2 + 2) - 2 * k * b21) - 0.375 * cn(5) * (3 * k^4 - 8 * k^2 + 8)) / (2 * lambda * (lambda * (1 + k^2) - 2 * k));
    s2 = (1.5 * cn(4) * (2 * a22 * (k^2 - 2) + a24 * (k^2 + 2) + 2 * k* b22 + 5 * d21) + 0.375 * cn(5) * (12 - k^2)) / (2 * lambda * (lambda * (1 + k^2) - 2 * k));
    
    a1 = -1.5 * cn(4) * (2 * a21 + a23 + 5 * d21) - 0.375 * cn(5) * (12 - k^2);
    a2 = +1.5 * cn(4) * (a24 - 2 * a22) + 1.125 * cn(5);
    l1 = a1 + 2 * lambda^2 * s1;
    l2 = a2 + 2 * lambda^2 * s2;

    deltan = won * n;

    % In-plane amplitude (related to Az by a non-linear analytical constraint)
    Ax = sqrt( (-del - l2 * Az^2) / l1 );

    % Phase space vector, third-order expansion
    x = a21 * Ax^2 + a22 * Az^2 - Ax * cos(tau) + (a23 * Ax^2 - a24 * Az^2) * cos(2*tau) + (a31 * Ax^3 - a32 * Ax * Az^2) * cos(3*tau);
    y = k * Ax * sin(tau) + (b21 * Ax^2 - b22 * Az^2) * sin(2*tau) + (b31 * Ax^3 - b32 * Ax * Az^2) * sin(3*tau);
    z = deltan* Az * cos(tau) + deltan * d21 * Ax * Az * (cos(2*tau) - 3) + deltan * (d32 * Az * Ax^2 - d31 * Az^3) * cos(3*tau);
    dx = +lambda * Ax * sin(tau) - 2 * lambda * (a23 * Ax^2 - a24 * Az^2) * sin(2*tau) - 3 * lambda * (a31 * Ax^3 - a32 * Ax * Az^2) *sin(3*tau);
    dy = +lambda * ( k * Ax * cos(tau) + 2 * (b21 * Ax^2 - b22 * Az^2) * cos(2*tau) + 3 * (b31 * Ax^3 - b32 * Ax * Az^2) * cos(3*tau) );
    dz = -lambda * deltan * Az * sin(tau) - 2 * lambda * deltan * d21 * Ax * Az * sin(2*tau) - 3 * lambda * deltan * (d32 * Az * Ax^2 - d31 * Az^3) * sin(3*tau);

    % Position vector
    r0 = [x; -y; z];           % Position vector
    r0 = gamma * r0;               % Re-scaled position vector

    % Velocity vector
    v0 = gamma * [dx; dy; dz];                            
    
    % Output
    seed = [r0; v0];
end