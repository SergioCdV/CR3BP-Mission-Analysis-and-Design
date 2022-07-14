%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 19/05/22

%% Plots %%
% Function to display the plots of the optimization

% Inputs: - structure system, containing the physical information of the
%           2BP of interest
%         - scalar tf, the final time of flight 
%         - vector tau, the final temporal grid 
%         - array C, the final state evolution matrix
%         - array u, a 3xm matrix with the control input evolution
%         - scalar T, the maximum allowed acceleration
%         - vector initial_coe, the initial orbital elements 
%         - vector final_coe, the final orbital elements 
%         - structure setup, containing the setup of the figures

function plots(system, tf, tau, C, u, T, initial_coe, final_coe, setup)
    % Set up 
    animations = setup.animations;

    % Constants 
    mu = system.mu;         % Gravitational parameter of the system 
    r0 = system.distance;   % Fundamental distance unit of the system
    t0 = system.time;       % Fundamental time unit of the system
    time = tf*tau; 

    % Setting up for gif image
    imsize = 350;
    im = cell(1,size(C,2));
    xpos = 10; ypos = 150;

    % Final state vector
    time_days = time*t0/86400;
    
    % Final spacecraft trajectory in Cartesian coordinates
    [S] = cylindrical2cartesian(C(1:3,:),true);
    x = S(1,:);
    y = S(2,:); 
    z = S(3,:);
    
    % Earth's orbit
    thetaE = linspace(0, 2*pi, size(C,2));
    
    s = coe2state(mu, initial_coe);
    initial = cylindrical2cartesian(s, false).';
    
    s = coe2state(mu, final_coe);
    final = cylindrical2cartesian(s, false).';
    
    s = zeros(6,length(thetaE));
    for i = 1:length(thetaE)
        s(:,i) = coe2state(mu, [initial_coe(1:end-1) initial(2)+thetaE(i)]);
    end
    xE = s(1,:);
    yE = s(2,:);
    zE = s(3,:);
    
    % Mars's orbit
    for i = 1:length(thetaE)
        s(:,i) = coe2state(mu, [final_coe(1:end-1) final(2)+thetaE(i)]);
    end
    xM = s(1,:);
    yM = s(2,:);
    zM = s(3,:);

    % Orbit representation
    figure_orbits = figure;
    view(3)
    set(figure_orbits,'position',[xpos,ypos,1.2*imsize,imsize])
    hold on
    xlabel('$X$ coordinate')
    ylabel('$Y$ coordinate')
    zlabel('$Z$ coordinate')
    plot3(0,0,0,'*k');
    plot3(x(1),y(1),z(1),'*k');
    plot3(xE,yE,zE,'LineStyle','--','Color','r','LineWidth',0.3);   % Earth's orbit
    plot3(xM,yM,zM,'LineStyle','-.','Color','b','LineWidth',0.3);   % Mars' orbit
    hold on
    grid on; 

    % Animation
    if (animations == 1)
        for i = 1:m
            marker = plot(x(i),y(i),'k*');
            hold on
            title(strcat('time elapsed = ',{' '}, num2str(round(time_days(i))),' days'));
            frame = getframe(figure_orbits);
            im{i} = frame2im(frame);
            delete(marker)
            plot(x(i),y(i),'k.');
        end
        
        frame = getframe(figure_orbits);
        im{i} = frame2im(frame);
        writegif('orbit.gif',im,m,2/m);
    end
    
    legend('off')
    %title('Optimal transfer orbit')
    plot3(x,y,z,'k','LineWidth',1);
    plot3(x(end),y(end),z(end),'*k');
    grid on;

    % Propulsive acceleration plot
    figure_propulsion = figure;
    set(figure_propulsion,'position',[xpos + 1.2*imsize,ypos,1.2*imsize,imsize])
    %title('Spacecraft acceleration in time')
    hold on
    plot(tau, sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)*r0/t0^2, 'k','LineWidth',1)
    plot(tau, u*r0/t0^2, 'LineWidth', 0.3)
    yline(T*r0/t0^2, '--k')
    xlabel('Flight time')
    ylabel('$\mathbf{a}$')
    legend('$a$','$a_\rho$','$a_\theta$','$a_z$')
    grid on;
    xlim([0 1])

    figure 
    hold on
    plot(time, rad2deg(atan2(u(2,:),u(1,:)))); 
    hold off 
    grid on;
    xlabel('Time')
    ylabel('$\theta$')
    title('Thrust in-plane angle')

    figure 
    hold on
    plot(time, rad2deg(atan2(u(3,:),sqrt(u(1,:).^2+u(2,:).^2)))); 
    hold off 
    grid on;
    xlabel('Time')
    ylabel('$\phi$')
    title('Thrust out-of-plane angle')
end