close all; clear all; clc;
%% parameters
% unit: length(m), mass(kg)
H = 1.9; % height of railing
l1 = 1.60; % height
l2 = 0.3; % length of down part of arm
l3 = 0.5; % length of up part of arm
L = 0.3;
m = 50; % weight
g = 9.8; %m/s^2
I = 1/3*m*l1^2; % moment of inertia
d0 = (L^2 + (H-l1)^2)^(1/2); % origin length of spring

% ----- time interval -----
t_react = 0.25;
t_return = 0.5; %settling time
% -------------------------

%% calculate acceleration(const)
v0 = 50/3.6; % m/s
dis = 20; % distance(m)
a = -v0^2/(2*dis);
t_total = -v0/a;

%% constraint
theta_max = asin(H/(H^2+L^2)^(1/2)) - asin((H^2+L^2+l1^2-l2^2-l3^2-2*l2*l3)/(2*l1*(H^2+L^2)^(1/2)));

%% reaction time
k = 0; % srping constant
c = 0; % damper constant
% ----- state space representation -----
sys = @(t,X)[X(2);Torque(X(1),X(2),l1,l2,l3,H,L,m,a,g,k,c,d0)/I];
% --------------------------------------

[ts_react,xs] = ode45(sys,[0,t_react],[0;0]);
theta_react = xs(:,1);    w_react = xs(:,2);
theta_index = min(find(theta_react > theta_max));
theta_react(theta_react>theta_max) = theta_max;
if (theta_index)
    w_react(w_react > w_react(theta_index-1,1)) = w_react(theta_index-1,1);
end
% plot(ts,t_result);

% ----- plot -----
% subplot(1,2,2);
% h1 = animatedline('Color','b','LineWidth',1);
% axis([min(ts_react(:,1)) max(ts_react(:,1)) min(theta_react) max(theta_react)]);
% title('theta-t plot');
% xlabel('t');
% ylabel('theta');
% 
% for i=1:size(ts_react,1)
%     subplot(1,2,1);
%     Plot(theta_react(i),H,L,l1,l2,l3);
%     addpoints(h1,ts_react(i,1),theta_react(i));
%     drawnow;
% end
% ----------------

%% balanced k,c
syms x
[O1,A1,B1,C1] = CalCoordinate(l1,l2,l3,H,L,theta_react(end,1));
equ = cross(1/2*[(A1-O1);0],m*(-a)*[-1;0;0]) + cross(1/2*[(A1-O1);0],m*g*[0;-1;0]) + cross([(A1-O1);0],(SpringForce(x,A1,C1,d0)+DampingForce(0,l1*w_react(end,1)))*[(C1-A1);0]/norm(C1-A1)) == [0;0;0];
k1 = double(vpa(solve(equ,x),7));

damp_r1 = 1; % damping ratio
c1 = 2*m*damp_r1*(k1/m)^(1/2);
t_bal = t_total - t_react;
% ----- state space representation -----
sys = @(t,X)[X(2);Torque(X(1),X(2),l1,l2,l3,H,L,m,a,g,k1,c1,d0)/I];
% --------------------------------------

[ts1,xs] = ode45(sys,[0,t_bal],[theta_react(end,1);w_react(end,1)]);
theta1 = xs(:,1);    w1 = xs(:,2);
theta1_index = min(find(theta1 > theta_max));
theta1(theta1 > theta_max) = theta_max;
if (theta1_index)
    w1(w1 > w1(theta_index-1,1)) = w1(theta_index-1,1);
end

% plot(ts,xs(:,1));
% ----- plot -----
% subplot(1,2,2);
% h1 = animatedline('Color','b','LineWidth',1);
% axis([min(ts(:,1)) max(ts(:,1)) min(xs(:,1)) max(xs(:,1))]);
% title('theta-t plot');
% xlabel('t');
% ylabel('theta');
% 
% for i=1:size(ts,1)
%     subplot(1,2,1);
%     Plot(xs(i,1),H,L,l1,l2,l3);
%     addpoints(h1,ts(i,1),xs(i,1));
%     drawnow;
% end
% ----------------

%% return to the original position
% syms x
% [O2,A2,B2,C2] = CalCoordinate(l1,l2,l3,H,L,theta1(end,1));
% equ = cross(1/2*[(A2-O2);0],m*a*[-1;0;0]) + cross(1/2*[(A2-O2);0],m*g*[0;-1;0]) + cross([(A2-O2);0],(SpringForce(x,A2,C2,d0)+DampingForce(0,l1*w1(end,1)))*[(C2-A2);0]/norm(C2-A2)) == [0;0;0];
% k2 = double(vpa(solve(equ,x),7));

damp_r2 = 0.6; % damping ratio
c2 = 2*m*(4/t_return);
k2 = m*(4/(t_return*damp_r2))^2;

a2 = 0; % car is stopped
% ----- state space representation -----
sys = @(t,X)[X(2);Torque(X(1),X(2),l1,l2,l3,H,L,m,a2,g,k2,c2,d0)/I];
% --------------------------------------

[ts2,xs] = ode45(sys,[0,t_return],[theta1(end,1);w1(end,1)]);
theta2 = xs(:,1);    w2 = xs(:,2);
theta2_index = min(find(theta2 > theta_max));
theta2(theta2 > theta_max) = theta_max;
if (theta2_index)
    w2(w2 > w2(theta_index-1,1)) = w2(theta_index-1,1);
end

%% assemble three stage data
theta = [theta_react ; theta1 ; theta2];
w = [w_react ; w1 ; w2];
ts = [ts_react ; ts1+t_react ; ts2+t_react+t_bal];

% ----- acceleration -----
acc = [a*ones(length([ts_react;ts1]),1) ; zeros(length(ts2),1)];
% ------------------------

% ----- velocity -----
vel = [v0.*ones(length([ts_react;ts1]),1) ; (v0+a*(t_react+t_bal))*ones(length(ts2),1)] + acc.*ts;
% vel = [v0*ones(length([ts_react;ts1]),1) + a.*[ts_react;ts1+t_react] ; zeros(length(ts2),1)];
% --------------------

%% static plot
figure(1)
subplot(4,1,1);
axis([min(ts(:,1)) max(ts(:,1)) min(vel) max(vel)]);
plot(ts(:,1),vel(:,1),'b','LineWidth',1);
axis tight;
title('v-t plot');
xlabel('t');
ylabel('v');

subplot(4,1,2);
axis([min(ts(:,1)) max(ts(:,1)) min(acc)-1 max(acc)]);
plot(ts(:,1),acc(:,1),'b','LineWidth',1);
axis tight;
title('acc-t plot');
xlabel('t');
ylabel('acc');

subplot(4,1,3);
axis([min(ts(:,1)) max(ts(:,1)) min(theta) max(theta)]);
plot(ts(:,1),theta(:,1),'b','LineWidth',1);
axis tight;
title('theta-t plot');
xlabel('t');
ylabel('theta');

subplot(4,1,4);
axis([min(ts(:,1)) max(ts(:,1)) min(w(:,1)) max(w(:,1))]);
plot(ts(:,1),w(:,1),'b','LineWidth',1);
axis tight;
title('w-t plot');
xlabel('t');
ylabel('w');

%% animation plot
% figure('units','normalized','outerposition',[0 0 1 1])
% filename = 'motion.gif';
% 
% subplot(4,2,2);
% h1 = animatedline('Color','b','LineWidth',1);
% axis([min(ts(:,1)) max(ts(:,1)) min(vel) max(vel)]);
% title('v-t plot');
% xlabel('t');
% ylabel('v');
% 
% subplot(4,2,4);
% h2 = animatedline('Color','b','LineWidth',1);
% axis([min(ts(:,1)) max(ts(:,1)) min(acc)-1 max(acc)]);
% title('acc-t plot');
% xlabel('t');
% ylabel('acc');
% 
% subplot(4,2,6);
% h3 = animatedline('Color','b','LineWidth',1);
% axis([min(ts(:,1)) max(ts(:,1)) min(theta) max(theta)]);
% title('theta-t plot');
% xlabel('t');
% ylabel('theta');
% 
% subplot(4,2,8);
% h4 = animatedline('Color','b','LineWidth',1);
% axis([min(ts(:,1)) max(ts(:,1)) min(w(:,1)) max(w(:,1))]);
% title('w-t plot');
% xlabel('t');
% ylabel('w');
% 
% del = 0.01; % time between animation frames
% for i=1:size(ts,1)
%     subplot(4,2,[1 3 5 7]);
%     Plot(theta(i,1),H,L,l1,l2,l3);
%     
%     addpoints(h1,ts(i,1),vel(i,1));
%     addpoints(h2,ts(i,1),acc(i,1));
%     addpoints(h3,ts(i,1),theta(i,1));
%     addpoints(h4,ts(i,1),w(i,1));
%     drawnow;
%     
%     % ----- save to gif -----
%     frame = getframe(2);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1;
%         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
%     end
% end

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
    xlim([-0.6 L+0.2]);
    ylim([0 H]);
end