close all; clear all; clc;
%% parameters
% unit: length(m), mass(kg)
H = 0.17; % height of railing
l1 = 0.141985; % height
l2 = 0.025; % length of down part of arm
l3 = 0.033; % length of up part of arm
L = 0.025;
m = 0.0848; % weight
g = 9.8; %m/s^2
% I = 0.000657;
I = 1/3*m*l1^2; % moment of inertia
d0 = (L^2 + (H-l1)^2)^(1/2); % origin length of spring

%% calculate acceleration(const)
v0 = 0.5; % m/s
% dis = 20; % distance(m)
% a = [-1 -2 -2.5 -3 -5];
a = -2;
% t_total = -v0/a;

%% constraint
theta_max = asin(H/(H^2+L^2)^(1/2)) - asin((H^2+L^2+l1^2-l2^2-l3^2-2*l2*l3)/(2*l1*(H^2+L^2)^(1/2)));
theta_min = -theta_max;

%% compare different acceleration
% damp_ratio = 0.16;
% k = 0.1237;
% c = 2*damp_ratio*I*(k/I)^(1/2);
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% for i = 1 : length(a)
%     ap = a(i);
%     
%     % ----- state space representation -----
%     sys = @(t,X)[X(2);Torque(X(1),X(2),l1,l2,l3,H,L,m,acc(ap,t),g,k,c,d0)/I];
%     % --------------------------------------
%     
%     [ts,xs] = ode45(sys,[0,2.2],[0;0]);
%     
%     theta = xs(:,1);
%     theta(theta>theta_max) = theta_max;
%     theta(theta<theta_min) = theta_min;
%     w = xs(:,2);    
%     
%     subplot(1,2,1);
%     plot(ts(:,1),theta(:,1),'LineWidth',1,'DisplayName',['a=',num2str(-ap)]); hold on;
%     axis([0 max(ts(:,1)) min(theta(:,1)) max(theta(:,1))]);
%     title('theta-t plot');
%     xlabel('t(s)');
%     ylabel('theta(rad)');
%     legend
%     
%     subplot(1,2,2);
%     plot(ts(:,1),w(:,1),'LineWidth',1,'DisplayName',['a=',num2str(-ap)]); hold on;
%     axis([0 max(ts(:,1)) min(w(:,1)) max(w(:,1))]);
%     title('w-t plot');
%     xlabel('t(s)');
%     ylabel('w(rad/s)');
%     legend
% end

%% experience
damp_ratio = 0.16;
k = 0.1237;
c = 2*damp_ratio*I*(k/I)^(1/2);
% c = 2*m*damp_ratio*(k/m)^(1/2);
% ----- state space representation -----
sys = @(t,X)[X(2);Torque(X(1),X(2),l1,l2,l3,H,L,m,acc(a,t),g,k,c,d0)/I];
% --------------------------------------

[ts,xs] = ode45(sys,[0,2.2],[0;0]);

% plot(ts(:,1),xs(:,1));

theta = xs(:,1);
theta(theta>theta_max) = theta_max;
theta(theta<theta_min) = theta_min;

w = xs(:,2);

%% animation plot
% ----- initialize video -----
% myVideo = VideoWriter(['torsionspring_',num2str(-a)]); %open video file
% myVideo.FrameRate = 24*4;
% open(myVideo);
% ----------------------------

figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,2,2);
h1 = animatedline('Color','b','LineWidth',1);
axis([0 max(ts(:,1)) min(theta(:,1)) max(theta(:,1))]);
title('theta-t plot');
xlabel('t');
ylabel('theta');

subplot(2,2,4);
h2 = animatedline('Color','b','LineWidth',1);
axis([0 max(ts(:,1)) min(w(:,1)) max(w(:,1))]);
title('w-t plot');
xlabel('t');
ylabel('w');

del = 0.01; % time between animation frames
for i=1:size(ts,1)
    subplot(2,2,[1,3]);
    Plot(theta(i,1),H,L,l1,l2,l3);
    
    addpoints(h1,ts(i,1),theta(i,1));
    addpoints(h2,ts(i,1),w(i,1));
    drawnow;
    
    % ----- save to video -----
%     frame = getframe(1);
%     writeVideo(myVideo,frame);
    
    % ----- save to gif -----
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    filename = ['torsionspring_',num2str(-a),'.gif'];
    if i == 1;
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
    end
end

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

function torque = Torque(theta,w,l1,l2,l3,H,L,m,a,g,kt,ct,d0)
    [O,A,B,C] = CalCoordinate(l1,l2,l3,H,L,theta);
    torque = cross(1/2*[(A-O);0],m*(-a)*[-1;0;0]) + cross(1/2*[(A-O);0],m*g*[0;-1;0]) + (kt*theta + ct*w)*[0;0;-1];
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

function ac = acc(a,t)
    if t < 0.1
        ac = a;
    else
        ac = 0;
    end
end