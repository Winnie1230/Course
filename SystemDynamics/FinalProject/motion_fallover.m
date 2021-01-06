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

% ----- friction coefficient -----
umax = 0.6; % maximum friction coefficient
% --------------------------------

%% calculate acceleration(const)
v0 = 60/3.6; % m/s
dis = 20; % distance(m)
a = -v0^2/(2*dis);
t_total = -v0/a;

%% constraint
theta_max = asin(H/(H^2+L^2)^(1/2)) - asin((H^2+L^2+l1^2-l2^2-l3^2-2*l2*l3)/(2*l1*(H^2+L^2)^(1/2)));

%% reaction time
k = 0; % spring constant
c = 0; % damper constant
% ----- state space representation -----
sys = @(t,X)[X(2);Torque(X(1),X(2),l1,l2,l3,H,L,m,a,g,k,c,d0)/I];
% --------------------------------------
figure('units','normalized','outerposition',[0 0 1 1])

[ts_react,theta_react,w_react] = React(sys,t_react,theta_max);
for i = 1:length(theta_react)
    Plot(theta_react(i,1),H,L,l1,l2,l3);
    drawnow;
end

% check if fallover
u = -a/g;
if u > umax
    iter = 50; %T/dt = 50
    [ts_fall,theta_fall,delta_x] = FallOver(sys,a,umax,g,l1,[theta_react(end);w_react(end)],iter);
    O = [0;0];
    [O,A,B,C] = CalCoordinate(l1,l2,l3,H,L,theta_react(end,1));
    % ----- plot fallover -----
    for i = 1 : length(theta_fall)
        O = O + [-1;0]*delta_x(i,1);
        PlotFallover(theta_fall(i,1),delta_x(i,1),O,l1,l2,l3,L,H);
        drawnow;
    end
    % --------------------------
    
else
    % ----- balanced k,c -----
    t_bal = t_total - t_react;
    ic = [theta_react(end,1);w_react(end,1)];
    [theta1,w1] = balanceKC(m,g,a,I,theta_max,l1,l2,l3,H,L,d0,t_bal,ic);
    for i = 1:length(theta1)
        Plot(theta1(i,1),H,L,l1,l2,l3);
        drawnow;
    end
    
    % ----- return to origin -----
    damp_r2 = 0.6; % damping ratio
    ic = [theta1(end,1);w1(end,1)];
    [theta2,w2] = ReturnOrigin(l1,l2,l3,H,L,m,g,I,d0,theta_max,damp_r2,t_return,ic);
end


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

function force = Force(k,A,C,d0,c,v)
    force = Spring(k,A,C,d0) + DampingForce(c,v);
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
%     xlim([-0.6 L+0.2]);
    xlim([-l1-1.5 L+0.2]);
    ylim([0 H]);
end

function PlotFallover(theta,dx,O,l1,l2,l3,L,H)
    A = O + l1*[cos(theta) -sin(theta) ; sin(theta) cos(theta)]*[0;1] + [-1;0]*dx;
    vec = (O-A)/norm(O-A);
    B = A + l2*vec;
    C = A + (l2+l3)*vec;
    plot([O(1),A(1)],[O(2),A(2)],'k','LineWidth',2); hold on;
%     plot([A(1),B(1)],[A(2),B(2)],'k','LineWidth',2); hold on;
%     plot([B(1),C(1)],[B(2),C(2)],'k','LineWidth',2); hold on;
    
    plot(O(1),O(2),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold on;
    plot(A(1),A(2),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold off;
%     plot(B(1),B(2),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold on;
%     plot(C(1),C(2),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6]); hold off;
    axis equal;
    xlim([-l1-1.5 L+0.2]);
    ylim([0 H]);
end

function [ts_react,theta_react,w_react] = React(sys,t_react,theta_max)
    [ts_react,xs] = ode45(sys,[0,t_react],[0;0]);
    theta_react = xs(:,1);    w_react = xs(:,2);
    theta_index = min(find(theta_react > theta_max));
    theta_react(theta_react>theta_max) = theta_max;
    if (theta_index)
        w_react(w_react > w_react(theta_index-1,1)) = w_react(theta_index-1,1);
    end
end

function [ts_fall,theta_fall,delta_x] = FallOver(sys,a,umax,g,l1,ic,iter)
    dt = 0.01;
    T_fall = 0:dt:dt*iter;
    [ts_fall,xs_fall] = ode45(sys,T_fall,ic);
    theta_fall = xs_fall(:,1);
    w_fall = xs_fall(:,2);
    fall_index = min(find(theta_fall > pi/2));
    theta_fall(fall_index : end,1) = pi/2;
    w_fall(fall_index:end,1) = 0;
    ax = -a - umax*g;
    v_fall = 1/2*l1*1/2.*(w_fall(1:end-1,1)+w_fall(2:end,1)) + ax*dt;
    delta_x(2:length(ts_fall),1) = v_fall*dt + 1/2*ax*dt^2;
    delta_x(1,1) = 0;
end

function [theta1,w1] = balanceKC(m,g,a,I,theta_max,l1,l2,l3,H,L,d0,t_bal,ic)
    syms x
    % react_end: theta_react(end,1);
    [O1,A1,B1,C1] = CalCoordinate(l1,l2,l3,H,L,ic(1));
    equ = cross(1/2*[(A1-O1);0],m*(-a)*[-1;0;0]) + cross(1/2*[(A1-O1);0],m*g*[0;-1;0]) + cross([(A1-O1);0],(SpringForce(x,A1,C1,d0)+DampingForce(0,l1*ic(2)))*[(C1-A1);0]/norm(C1-A1)) == [0;0;0];
    k1 = double(vpa(solve(equ,x),7));

    damp_r1 = 1; % damping ratio
    c1 = 2*m*damp_r1*(k1/m)^(1/2);

    % ----- state space representation -----
    sys = @(t,X)[X(2);Torque(X(1),X(2),l1,l2,l3,H,L,m,a,g,k1,c1,d0)/I];
    % --------------------------------------

    [ts1,xs] = ode45(sys,[0,t_bal],ic);
    theta1 = xs(:,1);    w1 = xs(:,2);
    theta1_index = min(find(theta1 > theta_max));
    theta1(theta1 > theta_max) = theta_max;
    if (theta1_index)
        w1(w1 > w1(theta_index-1,1)) = w1(theta_index-1,1);
    end
end

function [theta2,w2] = ReturnOrigin(l1,l2,l3,H,L,m,g,I,d0,theta_max,damp_r2,t_return,ic)
    c2 = 2*m*(4/t_return);
    k2 = m*(4/(t_return*damp_r2))^2;

    a2 = 0; % car is stopped
    % ----- state space representation -----
    sys = @(t,X)[X(2);Torque(X(1),X(2),l1,l2,l3,H,L,m,a2,g,k2,c2,d0)/I];
    % --------------------------------------

    [ts2,xs] = ode45(sys,[0,t_return],ic);
    theta2 = xs(:,1);    w2 = xs(:,2);
    theta2_index = min(find(theta2 > theta_max));
    theta2(theta2 > theta_max) = theta_max;
    if (theta2_index)
        w2(w2 > w2(theta_index-1,1)) = w2(theta_index-1,1);
    end
end
