function PhasePortrait(m,c,k,p0)
    % assume f(t) = 0, draw x-x_dot phase portrait using lots of initial conditions
    % ----- state space representation -----
    f2 = @(t,X)[X(2);-k/m*X(1)-c/m*X(2)];
    % --------------------------------------
    figure('Name','SMD system_phase portrait','units','normalized','outerposition',[0 0 1 1])

    for i = 1:length(p0)
        [ts,ys] = ode45(f2,[0 20],[p0(i);(-1)^i*p0(i)]);
        plot(ys(:,1),ys(:,2));  hold on
        plot(ys(1,1),ys(1,2),'bo'); hold on % starting point
        plot(ys(end,1),ys(end,2),'kx','MarkerSize',15,'LineWidth',4); hold on % ending point
        grid on
        axis equal tight;
    end
    hold off;
    title({['phase portrait(assume f(t)=0)'];['SMD system : m=',num2str(m),', c=',num2str(c),', k=', num2str(k)]});
    xlabel('$x$','interpreter','latex')
    ylabel('$\dot{x}$','interpreter','latex')
end