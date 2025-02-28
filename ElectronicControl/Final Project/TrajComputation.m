function TrajComputation(PQ,Q,target,TQ,iter,map_pos,map_len)
% use electronic field attractive force and repulsive force to let the
% particle get to the target position.

% This method has some problems:
% The attractive force become larger as the particle get closer to the
% target, so the acceleration become larger as well.
% As the attractive force become larger, the repulsive force of the obstacle may not be
% large enough to let the particle get away from the obstacle, and the
% planned trajectory may be collide with the obstacle.

    %% parameter
    q = 1; % car electron
    m = 0.001;
    k = 10;

    % target = [0;5];
    % TQ = -40;

    % T = 1;
    dt = 0.001;
    % iter = 1;

    % filename = 'motion.gif';
    % del = 0.01; % time between animation frames

    traj = zeros(iter,2);
    %% Trajectory computation
    % initial condition
    P = [0;0];
    v = [0;0];
    traj(1,:) = P;

    % Q = [2 1.5];
    % PQ = [-2,2 ; 1,0.5];

    % h = figure('units','normalized','outerposition',[0 0 1 1]);

    for i = 1:iter
        Fe = -(k*Q*q)'./vecnorm(P'-PQ,2,2).^3.*(PQ-P');
        Ft = -(k*TQ*q/norm(P-target)^3)*(target-P);

        Plot(P,Q,PQ,Fe,target,TQ,Ft); hold on
        plot(traj(1:i,1),traj(1:i,2),'Color','b','LineWidth',1); hold on
        PlotMap(map_pos,map_len); hold off
        axis equal
        axis([-50 50 0 100])
        drawnow

        F_total = sum(Fe,1)+Ft';
        a = F_total/m
        % next timestep velocity and position
        v = v + a'*dt;
        P = P + v*dt + 1/2*a'*dt

        norm(P-target)
        if norm(P-target) < 1
            break;
        end
        traj(i+1,:) = P';
        % ----- save to gif -----
    %     frame = getframe(h);
    %     im = frame2im(frame);
    %     [imind,cm] = rgb2ind(im,256);
    %     if i == 1
    %         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
    %     else
    %         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
    %     end
    end
end

%% function
function Plot(P,Q,PQ,Fe,target,TQ,Ft)
    factor = 10;
    
    circle(target(1), target(2), abs(TQ), [168/255 255/255 168/255]); hold on
    plot(target(1),target(2),'go','MarkerSize',10,'MarkerFace',[0 143/255 0]); hold on

    for i = 1:length(Q)
        circle(PQ(i,1), PQ(i,2), Q(i), [1 218/255 199/255]); hold on
        plot(PQ(i,1),PQ(i,2),'ro','MarkerSize',10,'MarkerFace','r'); hold on
    end
    for i = 1:length(Q)
        PlotQuiver(P,Fe(i,:)*factor,[245/255 78/255 0/255]); hold on
        PlotQuiver(PQ(i,:),-Fe(i,:)*factor,[245/255 78/255 0/255]); hold on
    end
    PlotQuiver(P,Ft*factor, [0/255 189/255 0/255]); hold on
    PlotQuiver(target,-Ft*factor,[0/255 189/255 0/255]); hold on
    plot(P(1),P(2),'bo','MarkerSize',10,'LineWidth',4); hold on
    hold off
end


function PlotQuiver(p1,dp,c)
%     dp = p2 - p1
    quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'Color',c,'MaxHeadSize',3);
end

function circle(x,y,r,c)
    th = 0:pi/50:2*pi;
    x_circle = r * cos(th) + x;
    y_circle = r * sin(th) + y;
    fill(x_circle, y_circle, c); hold on
    plot(x_circle, y_circle,'Color', c); hold off
%     axis equal
end
