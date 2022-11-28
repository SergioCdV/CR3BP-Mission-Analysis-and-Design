%% Autonomous RVD and docking in the CR3BP %% 
% Sergio Cuevas del Valle % 
% 28/11/22 % 

%% Numerical error %% 
% This script provides a script to test different nonlinear models for relative motion in the CR3BP. 

% The relative motion of two spacecraft in a halo orbit around L1 in the
% Earth-Moon system is analyzed both from the direct integration of the
% equations of motion and the full motion relative motion % 

% Units are non-dimensional and solutions are expressed in the synodic reference frame as defined by Howell, 1984.

%% Set up %%
close all; 
clear; 
clc; 

% Set path to Advanpix
addpath('C:\Users\sergi\Downloads\Multiprecision Computing Toolbox\');

% Set up graphics 
set_graphics();

%% Constants and initial data %% 
% MP test 
% mp.Test();

%% High accuracy integration %%
% Selection of precision 
mp.Digits(34);

% Conversion to desired precision 
mu = mp('0.0121505');
% s0 = mp('[0.8241 0 0.0568 0 0.1673 0 0.0006 0 0.0100 0 0.0114 0]');
load NumericalError.mat s0
s0 = mp(s0);
tspan = mp('0'):mp('1e-3'):mp('pi');

% Integration of the quadruple-precision relative models
options = odeset('RelTol', mp('2.25e-14'), 'AbsTol', repmat(mp('1e-22'),1,12));
[~, S_c] = ode113(@(t,s)dyn(mu, t, s), tspan, s0, options);

%% Low-accuracy integration
% Selection of precision 
mp.Digits(15);

% Conversion to desired precision 
mu = mp('0.0121505');
% s0 = mp('[0.8241 0 0.0568 0 0.1673 0 0.0006 0 0.0100 0 0.0114 0]');
load NumericalError.mat s0
s0 = mp(s0);
tspan = mp('0'):mp('1e-3'):mp('pi');

% Integration of the double-precision relative models
options = odeset('RelTol', mp('2.25e-14'), 'AbsTol', repmat(mp('1e-22'),1,12));
[~, S] = ode113(@(t,s)dyn_encke(mu, t, s), tspan, s0, options);
[~, Sn] = ode113(@(t,s)dyn(mu, t, s), tspan, s0, options);

%% Error comparison
error = S_c(:,7:12)-S(:,7:12);                  % Error via the Encke method
error_n = S_c(:,7:12)-Sn(:,7:12);               % Error via the full nonlinear model

e(:,1) = sqrt(dot(error, error, 2));            % State error (L2 norm) via Encke's method
e(:,2) = sqrt(dot(error_n, error_n, 2));        % State error (L2 norm) via Newtonian formulation

%% Plotting and results %% 
figure 
hold on
plot(tspan, log(e(:,1)), 'b'); 
plot(tspan, log(e(:,2)), 'r');
hold off
grid on
xlabel('$t$'); 
ylabel('Absolute error $\log{e}$');
legend('Encke', 'Newton');

%% Auxiliary functions
function [ds] = dyn(mu, t, s)
    % Preallocation 
    ds = mp(zeros(12,1));

    % Constants 
    mup(1) = mp('1')-mu; 
    mup(2) = mu; 

    R(:,1) = [-mu; mp('0'); mp('0')];
    R(:,2) = [mp('1')-mu; mp('0'); mp('0')];

    % State variables
    rt = s(1:3);
    vt = s(4:6);
    rho = s(7:9); 
    drho = s(10:12);

    rr(:,1) = rt-R(:,1);
    rr(:,2) = rt-R(:,2);

    % Absolute target motion 
    ds(1:3) = vt; 
    ds(4:6) = [rt(1)+2*vt(2); rt(2)-2*vt(1); mp('0')]-(mup(1)/norm(rr(:,1))^3*rr(:,1))-(mup(2)/norm(rr(:,2))^3*rr(:,2));

    % Relative motion 
    ds(7:9) = drho; 
    F(:,1) = mup(1)*(rr(:,1)/norm(rr(:,1))^3-(rho+rr(:,1))/norm(rho+rr(:,1))^3);     % Gravitational force of the first primary
    F(:,2) = mup(2)*(rr(:,2)/norm(rr(:,2))^3-(rho+rr(:,2))/norm(rho+rr(:,2))^3);     % Gravitational force of the second primary

    ds(10:12) = [rho(1)+2*drho(2); rho(2)-2*drho(1); mp('0')]+sum(F,2);
end

function [ds] = dyn_encke(mu, t, s)
    % Preallocation 
    ds = mp(zeros(12,1));

    % Constants 
    mup(1) = mp('1')-mu; 
    mup(2) = mu; 

    R(:,1) = [-mu; mp('0'); mp('0')];
    R(:,2) = [mp('1')-mu; mp('0'); mp('0')];

    % State variables
    rt = s(1:3);
    vt = s(4:6);
    rho = s(7:9); 
    drho = s(10:12);

    rr(:,1) = rt-R(:,1);
    rr(:,2) = rt-R(:,2);

    % Absolute target motion 
    ds(1:3) = vt; 
    ds(4:6) = [rt(1)+2*vt(2); rt(2)-2*vt(1); mp('0')]-(mup(1)/norm(rr(:,1))^3*rr(:,1))-(mup(2)/norm(rr(:,2))^3*rr(:,2));

    % Relative motion 
    ds(7:9) = drho; 

    gamma = [2*drho(2)+rho(1); -2*drho(1)+rho(2); mp('0')];
    for i = 1:length(mup)
        q = -dot(2*rr(:,i)+rho,rho)/norm(rho+rr(:,i))^2;
        f = q*(3*(1+q)+q^2)/(1+(1+q)^(3/2));
        gamma = gamma - (mup(i)/norm(rr(:,i))^3)*(f*rr(:,i)+(1+f)*rho);
    end

    ds(10:12) = gamma;
end
