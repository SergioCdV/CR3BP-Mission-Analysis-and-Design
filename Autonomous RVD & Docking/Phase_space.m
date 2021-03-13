%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 13/03/21 % 

%% Phase space study around the leader %% 
 
% The eigenvalues of the full linear model around the leader are solved analytically % 

% Units are non-dimensional and solutions are expressed in the Lagrange
% points reference frame as defined by Howell, 1984.

%% Initial data
%Gravitational parameters 
mu1 = @(mu)1-mu;                 %Gravitational parameter of the first primary
mu2 = @(mu)mu;                   %Gravitational parameter of the second primary

%Position of the primaries  
R1 = @(mu)[-mu; 0; 0];           %Position of the first primary
R2 = @(mu)[1-mu; 0; 0];          %Position of the second primary 

%Target position 
r = @(x,y,z)[x; y; z];

%Relative position to the primaries 
r1 = @(mu,x,y,z)r-R1;            %Relative position vector of the target to the first primary
r2 = @(mu,x,y,z)r-R2;            %Relative position vector of the target to the second primary
ur1 = @(mu,x,y,z)r1/norm(r1);    %Relative unit position vector of the target to the first primary
ur2 = @(mu,x,y,z)r2/norm(r2);    %Relative unit position vector of the target to the second primary

%State matrix 
O = zeros(3,3); 
I = eye(3); 
H = @(mu,x,y,z)-((mu1)/(norm(r1)^3)+(mu2)/(norm(r2)^3))*eye(3)+3*((mu1)/(norm(r1)^3)*(ur1*ur1.')+(mu2)/(norm(r2)^3)*(ur2*ur2.'));
S = [0 2 0; 0 -2 0; 0 0 0];
A = @(mu,x,y,z)[O I; H S]; 

%Eigenspectrum
[W, e] = eig(@(mu,x,y,z)A);
