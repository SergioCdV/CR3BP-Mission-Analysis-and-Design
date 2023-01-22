%% Project: 
% Date: 29/01/22

%% Approximating general function by orthornomal Bernstein polynomials
% This script approximates general 3D curves by means of orthogonal Berstein
% poylnomials

clear
close all

%% General setup 
% Define number of steps desired
steps = 1e4;

%% Generate the function 
% Curve parametrization
t0 = -pi;                        % Initial parameter
tf = pi;                         % Final parameter
t = linspace(t0, tf, steps);    % Curve parameter

tvec = (t-t(1))/(t(end)-t(1));


% 3D curve definition 
C = [sin(4*t); cos(4*t); sin(3*t)];

%% Function approximation
% Order of curve
n = 20;

% Computation of the orthogonal Berstein polynomials 
B = OB_basis(tvec,n);

% Compute the control points
P = zeros(size(C,1),n+1);
for i = 1:n
    for j = 1:size(C,1)
        P(j,i) = trapz(t,B(i,:).*C(j,:))/(t(end)-t(1));   
    end
end

% Generate the approximation
B1 = OB_classic(P,tvec);

for j = 1:size(C,1)
    P(j,:) = C(j,:)*pinv(B);
end

% Generate the approximation
B2 = OB_classic(P,tvec);

% Plot results
figure
view(3)
hold on
plot3(C(1,:), C(2,:), C(3,:));
plot3(B1(1,:), B1(2,:), B1(3,:));
plot3(B2(1,:), B2(2,:), B2(3,:),'k');
hold off
grid on 
legend('Original curve', 'Bézier approximation', 'LS Bézier approximation', 'Location', 'northwest')
title('Function approximation by means of orthogonal Bernstein polynomials');
