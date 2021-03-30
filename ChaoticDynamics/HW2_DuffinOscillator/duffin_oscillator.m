clear all; close all; clc;
% ----- parameters -----
delta = 0.25;
beta = 1;
w = 1.4;
gamma = 0.3;
% ----------------------

% ----- state space representation -----
f = @(t,X)[X(2) ; beta*X(1)-X(1)^3-delta*X(2)+gamma*cos(w*t)];
% --------------------------------------
ic = [1;1];

%% phase portrait
figure('units','normalized','outerposition',[0 0 1 1])
truncate_phase = 1000;
[ts,ys] = ode45(f,[0 1000],ic);
plot(ys(truncate_phase:end,1),ys(truncate_phase:end,2),'b');
grid on
axis equal tight;
axis([-2 2 -2 2]);

title({['phase portrait'];['delta=',num2str(delta),', beta=',num2str(beta),', w=', num2str(w),', gamma=',num2str(gamma)]});
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')

%% poincare plot
figure('units','normalized','outerposition',[0 0 1 1])
n = 1;
dt = n*2*pi/w;
iter = 2000;
T = 0 : dt : iter*dt;
trun_poincare = 30; % truncate point of poincare plane
[ts,ys] = ode45(f,T,ic);

plot(ys(trun_poincare:end,1),ys(trun_poincare:end,2),'b*');
title({['Poincare Plot of Duffin''s Oscillator'];['delta=',num2str(beta),', beta=',num2str(beta),', w=', num2str(w),', gamma=',num2str(gamma),', n=',num2str(n)]});
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')
grid on
axis equal tight;
axis([-2 2 -2 2]);