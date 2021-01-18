%% Itinerary computation in the CR3BP %%
% Sergio Cuevas del Valle, 26/09/20

%% Set up conditions 
setup_path();
format long;
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-14); 

%% Energy level and the 3 body sytem
mu = 3.054248395728148e-06;
eps = -4;

%% Lagrange point position 
L = libration_points(mu); 
xe = L(1,1);

%% Family of periodic orbits
A = 1e-4;
bmu = mu*abs(xe-1+mu)^(-3)+(1-mu)*abs(xe+mu)^(-3);
nu = sqrt(-0.5*(bmu-2-sqrt(9*bmu^2-8*bmu)));
tau = -(nu^2+2*bmu+1)/(2*nu);
vy0 = -A*tau*nu;
p0 = xe-A;    
phi0 = eye(4);
dx = 1e-9;

figure(1)
hold on
for i = 1:1
    %Initial conditions (continuation)
    p0 = p0+dx*(i-1);
    x0(i,1,:) = [p0 0 0 vy0 phi0(:).'];

    %Differential correction scheme (simple targetting)
    saveSol(1) = 10^8;
    tol = 4e-5; 
    vTol = 1e-8; 
    itermax = 10;
    j = 2;
    go_on = true;

    while (go_on) && (j < itermax)
        %Dynamics integration 
        [sol, phi, elaps] = integration_loop(mu, shiftdim(x0(i,end,:)));
    
        %Save the solution 
        saveSol(j) = sol(end,3);
    
        %Tolerance reaching supervision 
        if (abs(sol(end,3)) < vTol)
            go_on = false;
        else
            %Differential correction
            dv = sol(end,3)/(phi(3,4)-(sol(end,3)/sol(end,4))*phi(2,4));
            x0(i,j,:) = x0(i,end,:); 
            x0(i,j,4) = x0(i,j,4)-dv;
            j = j+1;
        end
    end
    
    %Compute the full periodic orbit
    dt = 1e-4; 
    tspan = 0:dt:2*elaps;
    [~, r] = ode113(@(t,x)cr3bp_dynamics(mu, x, 1), tspan, x0(i,end,:), options);

    %Plot the orbit 
    plot(r(:,1), r(:,2));        %Compute the full periodic orbit
    dt = 1e-4; 
    tspan = 0:dt:2*elaps;
    [~, r] = ode113(@(t,x)cr3bp_dynamics(mu, x, 1), tspan, x0(i,end,:), options);

    %Plot the orbit 
    plot(r(:,1), r(:,2));
end
plot(xe,0,'ko');
plot(1-mu,0,'ro');
hold off
        
%% Associated invariant manifolds
epsilon = 1e-8;
%Compute the monodromy matrix
k = 0;
for i = 1:4
    for j = 1:4
        mPhi(i,j) = r(end,5+k);
        k = k+1;
    end
end 
[W, eigs] = eig(mPhi);
x0 = x0(end,:);
T = 2*elaps;
tf = 3*T;
tspan = linspace(0, tf, 2000);
points = 50; 

%Select HOI points from the unstable/stable manifold to fiber it 
h = round((size(r,1)-1)/points); 
for i = 1:points
    haloX0(i,:) = r((i-1)*h+1,1:4); 
    haloT(i) = (i-1)*h+1;
end

%Select the unstable/stable vector (test case: unstable) 
Vs = W(:,4)/norm(W(:,4));
M0 = [haloX0(1,:)+epsilon*Vs.'];

%Generate the initial conditions along the fiber
for i = 2:points
    l = 0;
    for j = 1:4
        for k = 1:4
            Phi(j,k) = r(haloT(points), 5+l); 
            l = l+1;
        end
    end
    Vs = Phi*W(:,1);
    Vs = Vs/norm(Vs);
    M0(i,:) = [haloX0(i,:)+epsilon*Vs.'];
end

%Compute the tube and plot results
figure(2)
hold on
for i = 1:points 
    [t, M] = ode113(@(t,x)cr3bp_equations(mu, x, 1), tspan, M0(i,:), options);
    plot(M(:,1), M(:,2), 'r');
end
plot(r(:,1), r(:,2));
plot(xe,0,'ro');
plot(1-mu,0,'ko')
hold off
title('Unstable manifold');

%% Auxiliary functions 
function [x, phi, elaps] = integration_loop(mu, x0)
    %Set integration tolerances 
    RelTol = 1e-10;
    AbsTol = 1e-10;
    options = odeset('RelTol', RelTol, 'AbsTol', AbsTol, 'Events', @(t,x)x_crossing(t,x));

    %Set integration parameters
    x0 = x0.';
    yTol = 1e-11;
    elaps = 0;
    dt = 1e-3;  
    go_on = true; 
    crossing = false;
    i = 1; 

    while (go_on)
        %Update the tiem step to refine the results
        Dt = [0 dt];
        
        %Integration 
        [~,x_aux,~,~,eventIndex] = ode113(@(t,x)cr3bp_dynamics(mu, x, 0), Dt, x0(end,:), options);
        
        %Update initial conditions
        x0(i+1,:) = x_aux(end,:);

        %Tolerance reaching supervision
        if (eventIndex == 1)
            go_on = false;
            crossing = true; 
            dt = 1e-8;             
        elseif ((crossing) && (abs(x0(end,2)) < yTol))
            go_on = false;
        end
        
        %Elapsed time
        elaps = elaps+dt;
        
        %Next iteration counter update
        i = i+1;
    end
    
    %Final state evolution 
    x = x0;
    
    %Compute the state transition matrix 
    k = 0;
    for i = 1:4
        for j = 1:4
            phi(i,j) = x(end,5+k);
            k = k+1;
        end
    end
end

function [dx] = cr3bp_dynamics(mu, p, direction)
    %State vector
    x = p(1); 
    y = p(2);  
    v = p(3:4); 
    
    %Relevant parameters
    mu1 = 1-mu;
    mu2 = mu;
    r1 = [(x+mu2); y];
    r2 = [(x-mu1); y];
    
    %CR3BP Dynamics 
    if (direction == -1)
        inert = [x-2*v(2); y+2*v(1)];
    else
        inert = [x+2*v(2); y-2*v(1)];
    end
    
    dx = [v; inert-mu1/norm(r1)^3*r1-mu2/norm(r2)^3*r2];
    
    %Compute the state transition matrix evolution
    k = 0;
    for i = 1:4
        for j = 1:4
            phi0(i,j)= p(5+k);
            k = k+1;
        end
    end
    
    %Compute the variational equations
    Ux = -1+(mu1/norm(r1)^3)*(1-3*((x+mu2)/(norm(r1)))^2)+(mu2/norm(r2)^3)*(1-3*((x-mu1)/(norm(r2)))^2); 
    Uy = -1+(mu1/norm(r1)^3)*(1-3*(y/(norm(r1)))^2)+(mu2/norm(r2)^3)*(1-3*(y/(norm(r2)))^2);
    Uyx = -3*((mu1/norm(r1)^5)*y*(x+mu2)+(mu2/norm(r2)^5)*y*(x-mu1));
    Uxy = Uyx;
    Jacob = [0 0 1 0; 0 0 0 1; -Ux -Uxy 0 2; -Uyx -Uy -2 0];
    dphi = (Jacob*phi0).';
    
    %Update the whole state
    dx = [dx; dphi(:)];
    
    if (direction == -1)
        dx = -dx;
    end
end

function [dx] = cr3bp_equations(mu, p, direction)
    %State vector
    x = p(1); 
    y = p(2);  
    v = p(3:4); 
    
    %Relevant parameters
    mu1 = 1-mu;
    mu2 = mu;
    r1 = [(x+mu2); y];
    r2 = [(x-mu1); y];
    
    %CR3BP Dynamics 
    if (direction == -1)
        inert = [x-2*v(2); y+2*v(1)];
    else
        inert = [x+2*v(2); y-2*v(1)];
    end
    
    dx = [v; inert-mu1/norm(r1)^3*r1-mu2/norm(r2)^3*r2];
    
    if (direction == -1)
        dx = -dx; 
    end
end

function [position,isterminal,direction] = orbitEvents(t,x)   
    %Event defining
    position = x(2); 
    isterminal = 1; 
    direction = -1; 
end