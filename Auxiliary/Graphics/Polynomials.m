%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 24/02/21
% File: polynomials_graphics.m 
% Issue: 0 
% Validated: 

%% Polynomials graphics %% 
% This script provides a function to plot several Legendre and Chebyshev
% polynomials 

%% Main computation %%
set_graphics(); 

%Initial data
order = 5;                                          %Number of polynomials to plot
u = -1:1e-2:1;                                      %Normalized domain

%Polynomials evaluation
L = zeros(order, length(u));                        %Preallocation of the polynomial basis
T = zeros(order, length(u));                        %Preallocation of the polynomial basis

for i = 1:length(u)
    L(:,i) = legendre_polynomials(order, u(i));     %Computation of the polynomials
    T(:,i) = chebyshev('first', order, u(i));       %Computation of the polynomials
end

%% Plots and results %% 
figure(1) 
hold on 
for i = 1:order
    plot(u, L(i,:));
end
hold off
grid on; 
xlabel('Normalized domain $u$');
legend('$P_0(u)$', '$P_1(u)$', '$P_2(u)$', '$P_3(u)$', '$P_4(u)$')
ylabel('Legendre polynomial $P_n(u)$');

figure(2) 
hold on 
for i = 1:order
    plot(u, T(i,:));
end
hold off
grid on; 
xlabel('Normalized domain $u$');
% legend('$T_0(u)$', '$T_1(u)$', '$T_2(u)$', '$T_3(u)$', '$T_4(u)$');
ylabel('Chebyshev polynomial $T_n(u)$');

