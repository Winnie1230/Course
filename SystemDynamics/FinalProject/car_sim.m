close all; clear all; clc;
%% parameters
% unit: length(m), mass(kg)
H = 0.17; % height of railing
l1 = 0.141985; % height
l2 = 0.022; % length of down part of arm
l3 = 0.033; % length of up part of arm
L = 0.025;
m = 0.0848; % weight
g = 9.8; %m/s^2
I = 0.000657;
% I = 1/3*m*l1^2; % moment of inertia
d0 = (L^2 + (H-l1)^2)^(1/2); % origin length of spring

u = 0.6; % friction coefficient

damp_ratio = 0.26;
k = 87;
c = 2*I*damp_ratio*(l1^2*k/I)^(1/2)/l1^2;

t_last = 0.1; % time that accelerate affect
t_sim = 1; % time of total simulation time

%% calculate acceleration(const)
v0 = 0.5; % m/s
% dis = 20; % distance(m)
% a = [-1 -2 -2.1];
% a = [-1 -2 -3 -3.5 -3.9 -4];
a = -3.5;
% t_total = -v0/a;

%% constraint
theta_max = asin(H/(H^2+L^2)^(1/2)) - asin((H^2+L^2+l1^2-l2^2-l3^2-2*l2*l3)/(2*l1*(H^2+L^2)^(1/2)));
theta_min = -theta_max;

%% compare different acceleration
plot_num = 3;

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1 : length(a)
    ap = a(i);
    % ----- state space representation -----
    sys = @(t,X)[X(2);Torque(X(1),X(2),l1,l2,l3,H,L,m,acc(ap,t,t_last),g,k,c,d0)/I];
    % --------------------------------------

    [ts,xs] = ode45(sys,[0,t_sim],[0;0]);

    theta = xs(:,1);
    theta(theta>theta_max) = theta_max;
    theta(theta<theta_min) = theta_min;
    w = xs(:,2);
    
    % fk, fs_max
    fk = zeros(length(ts(:,1)),1);
%     fspring = zeros(length(ts(:,1)),1);
%     fdamp = zeros(length(ts(:,1)),1);
%     fa = zeros(length(ts(:,1)),1);
    fs = zeros(length(ts(:,1)),1);
    fs_max = zeros(length(ts(:,1)),1);
    flag = 0;
    for i = 1:length(ts(:,1))
        [O,A,B,C] = CalCoordinate(l1,l2,l3,H,L,theta(i,1));
%         fs = cross([(A-O);0],SpringForce(k,A,C,d0)*[(C-A);0]/norm(C-A));
%         fspring(i,1) = fs(3);
%         fd = cross([(A-O);0],(DampingForce(c,l1*w(i,1)))*[(C-A);0]/norm(C-A));
%         fdamp(i,1) = fd(3);
%         faa = cross(1/2*[(A-O);0],m*(-a)*[-1;0;0]);
%         fa(i,1) = faa(3);
        fk(i,1) = SpringForce(k,A,C,d0)+DampingForce(c,l1*w(i,1));
        fs(i,1) = m*(-ap) - fk(i,1)*cos(theta(i,1)); % friction
        fs_max(i,1) = sign(fs(i,1))*u*(m*g-fk(i,1)*sin(theta(i,1))); % friction maximum
        if abs(fs(i,1)) > abs(fs_max(i,1))
            fs(i,1) = fs_max(i,1);
            theta(i+1:end,1) = 0;
            w(i+1:end,1) = 0;

            flag = 1;
            break;
        end
    end
    
    subplot(plot_num,1,1);
    plot(ts(:,1),theta(:,1),'LineWidth',1,'DisplayName',['a=',num2str(-ap)]); hold on;
    title({['damping ratio=' num2str(damp_ratio),', k=',num2str(k),', c=',num2str(c)];['theta-t plot']});
    xlabel('t(s)');
    ylabel('theta(rad)');
    legend
    grid on;
    
    subplot(plot_num,1,2);
    plot(ts(:,1),w(:,1),'LineWidth',1,'DisplayName',['a=',num2str(-ap)]); hold on;
    title('w-t plot');
    xlabel('t(s)');
    ylabel('w(rad/s)');
    legend
    grid on;
    
    subplot(plot_num,1,3);
    plot(ts(:,1),fs(:,1),'LineWidth',1,'DisplayName',['a=',num2str(-ap)]); hold on;
    title('friction-t plot');
    xlabel('t(s)');
    ylabel('friction(N)');
    legend
    grid on;
    
%     subplot(plot_num,1,4);
%     plot(ts(:,1),fspring(:,1),'LineWidth',1,'DisplayName',['a=',num2str(-ap)]); hold on;
%     title('springforce-t plot');
%     xlabel('t(s)');
%     ylabel('spring_force(N)');
%     legend
%     grid on;
%     
%     subplot(plot_num,1,5);
%     plot(ts(:,1),fdamp(:,1),'LineWidth',1,'DisplayName',['a=',num2str(-ap)]); hold on;
%     title('dampingforce-t plot');
%     xlabel('t(s)');
%     ylabel('dampingforce(N)');
%     legend
%     grid on;
%     
%     subplot(plot_num,1,6);
%     plot(ts(:,1),fa(:,1),'LineWidth',1,'DisplayName',['a=',num2str(-ap)]); hold on;
%     title('ma-t plot');
%     xlabel('t(s)');
%     ylabel('ma');
%     legend
%     grid on; 
    
    if flag == 1
        subplot(plot_num,1,1);
        t1 = text(ts(i,1),theta(i,1),'\leftarrow slip away');
        t1.Color = 'red';
        
        subplot(plot_num,1,2);
        t = text(ts(i,1),w(i,1),'\leftarrow slip away');
        t.Color = 'red';    
        
        subplot(plot_num,1,3);
        t = text(ts(i,1),fs(i,1),'\leftarrow slip away');
        t.Color = 'red';
    end
end

%% experience
% % ----- state space representation -----
% sys = @(t,X)[X(2);Torque(X(1),X(2),l1,l2,l3,H,L,m,acc(a,t,t_last),g,k,c,d0)/I];
% % --------------------------------------
% 
% [ts,xs] = ode45(sys,[0,t_sim],[0;0]);
% 
% theta = xs(:,1);
% theta(theta>theta_max) = theta_max;
% theta(theta<theta_min) = theta_min;
% 
% w = xs(:,2);

%% animation plot
% initialize video
% myVideo = VideoWriter('myVideoFile'); %open video file
% myVideo.FrameRate = 24*4;
% open(myVideo);

% figure('units','normalized','outerposition',[0 0 1 1])
% 
% subplot(2,2,2);
% h1 = animatedline('Color','b','LineWidth',1);
% axis([0 max(ts(:,1)) min(theta(:,1)) max(theta(:,1))]);
% title('theta-t plot');
% xlabel('t');
% ylabel('theta');
% 
% subplot(2,2,4);
% h2 = animatedline('Color','b','LineWidth',1);
% axis([0 max(ts(:,1)) min(w(:,1)) max(w(:,1))]);
% title('w-t plot');
% xlabel('t');
% ylabel('w');
% 
% del = 0.01; % time between animation frames
% for i=1:size(ts,1)
%     subplot(2,2,[1,3]);
%     Plot(theta(i,1),H,L,l1,l2,l3);
% %     title({['car simulation -- a=',num2str(-a)]});
%     
%     addpoints(h1,ts(i,1),theta(i,1));
%     addpoints(h2,ts(i,1),w(i,1));
%     drawnow;
%     
%     % ----- save to video -----
% %     frame = getframe(1);
% %     writeVideo(myVideo,frame);
%     
%     % ----- save to gif -----
% %     frame = getframe(1);
% %     im = frame2im(frame);
% %     [imind,cm] = rgb2ind(im,256);
% %     filename = ['spring_',num2str(-a),'.gif'];
% %     if i == 1;
% %         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
% %     else
% %         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
% %     end
% end

% close(myVideo);

%% function
function [O,A,B,C] = CalCoordinate(l1,l2,l3,H,L,theta)
    O = [0;0];
    A = l1*[cos(theta) -sin(theta) ; sin(theta) cos(theta)]*[0;1];
    C = [L;H];
    d = norm(A-C); % spring length 
    
    beta = atan2(H-l1*cos(theta),L+l1*sin(theta));
    theta1 = acos((d^2+l2^2-l3^2)/(2*d*l2));
    alpha = beta - theta1;

    B = A + l2*[cos(alpha) -sin(alpha) ; sin(alpha) cos(alpha)]*[1;0];
end

function torque = Torque(theta,w,l1,l2,l3,H,L,m,a,g,k,c,d0)
    [O,A,B,C] = CalCoordinate(l1,l2,l3,H,L,theta);
    torque = cross(1/2*[(A-O);0],m*(-a)*[-1;0;0]) + cross(1/2*[(A-O);0],m*g*[0;-1;0]) + cross([(A-O);0],(SpringForce(k,A,C,d0)+DampingForce(c,l1*w))*[(C-A);0]/norm(C-A));
    torque = torque(3);
end

function fs = SpringForce(k,A,C,d0)
    fs = k*(norm(C-A)-d0);
end

function fd = DampingForce(c,v)
    fd = c*v;
end

function Plot(theta,H,L,l1,l2,l3)
    [O,A,B,C] = CalCoordinate(l1,l2,l3,H,L,theta);
    plot([O(1),A(1)],[O(2),A(2)],'k','LineWidth',2); hold on;
    plot([A(1),B(1)],[A(2),B(2)],'k','LineWidth',2); hold on;
    plot([B(1),C(1)],[B(2),C(2)],'k','LineWidth',2); hold on;
    
    plot(O(1),O(2),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold on;
    plot(A(1),A(2),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold on;
    plot(B(1),B(2),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold on;
    plot(C(1),C(2),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold off;
    
    axis equal;
    xlim([-0.05 L+0.05]);
    ylim([0 H]);
end

function ac = acc(a,t,t_last)
    if t < t_last
        ac = a;
    else
        ac = 0;
    end
end