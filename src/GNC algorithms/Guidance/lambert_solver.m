%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 12/10/21
% File: lambert.m 
% Issue: 0 
% Validated: 

%% Lambert solver %%
% This function contains a basic Lambert solver as given in Curtis, Orbital Mechanics for Engineering Students.

% Inputs: - scalar mu, the reduced gravitational parameter of the system
%         - vector r1, the initial position vector
%         - vector r2, the final position vector
%         - tf, the required time of flight

% Outputs: - vector v1 and v2, the terminal velocities at r1 and r2
%          - vector elements, the classical orbital elements of the
%            computed orbit

% Methods: 

function [v1, v2, path] = lambert_solver(mu, r1, r2)
    %Position vector norms 
    R1 = norm(r1);          %Norm of the initial position vector
    R2 = norm(r2);          %Norm of the final position vector
    r_1m = R1; 
    r_2m = R2; 

    %Anomaly desambiguation 
    z = cross(r1, r2);                  %Normal vector to the orbital plane
    cdtheta = dot(r1,r2)/(R1*R2);       %Cosine of the anomaly

    if (z(3) < 0)
        dtheta = 2*pi-acos(cdtheta);    %Change in the anomaly
    else
        dtheta = acos(cdtheta);         %Change in the anomaly
    end

    %Selection of the shortest path 
    if (dtheta > pi)
        tm = -1; 
        path = 'Long path';
    else
        tm = 1; 
        path = 'Short path';
    end

    %Sine of the anomaly change 
    sdtheta = tm*sqrt(1-cdtheta^2);

    %Main computation usign the f and g coefficients and a Newton method 
    iter = 1;               %Initial iteration
    GoOn = true;            %Convergence boolean
    maxIter = 1e3;          %Maximum number of iterations
    tol = 1e-3;             %Convergence tolerance 
    k = 5;                  %Method coefficient 
    
    %Preallocation of the angular momentum norm 
    h = zeros(1,maxIter);   %Angular momentum norm
    h(1) = norm(z);         %Initial iteration

    while ((iter < maxIter) && (GoOn)) 
        %Coefficients of the functions 
        f_a = 1-(1-cdtheta)*(mu*r_2m/h(iter)^2); 
        df_a = (1-cdtheta)*(2*mu*r_2m/h(iter)^3);
        f_b = ((1-cdtheta)/sdtheta)*(mu/h(iter))*((1-cdtheta)*(mu/h(iter)^2)-(1/r_1m)-(1/r_2m));
        df_b = ((1-cdtheta)/sdtheta)*((-mu/h(iter)^2)*((1-cdtheta)*(mu/h(iter)^2)-(1/r_1m)-(1/r_2m))+(mu/h(iter))*((1-cdtheta)*(-2*mu/h(iter)^3)));
        f_c = (r_1m*r_2m*sdtheta)/h(iter); 
        df_c = (-r_1m*r_2m*sdtheta)/h(iter)^2;
        f_d = 1-(1-cdtheta)*(mu*r_1m/h(iter)^2); 
        df_d = (1-cdtheta)*(2*mu*r_1m/h(iter)^3);

        ddf_a = (1-cdtheta)*(-6*mu*r_2m/h(iter)^4);
        ddf_b1 = ((1-cdtheta)/sdtheta)*((2*mu/h(iter)^3)*((1-cdtheta)*(mu/h(iter)^2)-(1/r_1m)-(1/r_2m))+(-mu/h(iter)^2)*((1-cdtheta)*(-2*mu/h(iter)^3))); 
        ddf_b2 = ((1-cdtheta)/sdtheta)*((-mu/h(iter)^2)*((1-cdtheta)*(-2*mu/h(iter)^3))+(mu/h(iter))*((1-cdtheta)*(6*mu/h(iter)^4))); 
        ddf_b = ddf_b1+ddf_b2;
        ddf_c = (2*r_1m*r_2m*sdtheta)/h(iter)^3;
        ddf_d = (1-cdtheta)*(-6*mu*r_1m/h(iter)^4);

        %Derivative of the problem function 
        f = f_a*f_d-f_b*f_c-1; 
        df = (df_a*f_d+f_a*df_d)-(df_b*f_c+f_b*df_c);
        ddf = (ddf_a*f_d+2*df_a*df_d+f_a*ddf_d)-(ddf_b*f_c+2*df_b*df_c+f_b*ddf_c);

        %Newton method iteration 
        ds = -k*f/(df+sqrt((k-1)^2*df^2-k*(k-1)*f*ddf));
        h(iter+1) = h(iter) + ds; 

        %Convergence analysis 
        if (abs(ds) < tol)
            GoOn = false; 
        else
            iter = iter+1; 
        end
    end

    %Final orbit computation 
    h = h(iter); 
    f = 1-(1-cdtheta)*(mu*r_2m/h(end)^2);      
    df = ((1-cdtheta)/(sdtheta))*(mu/h)*((mu/h^2)*(1-cdtheta)-1/r_1m-1/r_2m);
    g = (r_1m*r_2m*sdtheta)/h;
    dg = 1-(1-cdtheta)*(mu*r_1m/h^2);

    %Final velocity vectors
    v1 = (1/g)*(r2-f*r1);       %Initial velocity vector
    v2 = (1/g)*(dg*r2-r1);      %Final velocity vector
end