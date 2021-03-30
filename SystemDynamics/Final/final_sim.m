clear all; close all; clc;
% ----- parameters -----
global mb rg lc l mp k c g
mb = 1;
rg = 1;
lc = 1;
l = 2;
mp = 1;
k = 1;
c = 1;
g = 9.8;

ic = [1;0;0;0]; % initial condition
T = @(t)(0);
% T = @(t)(2*sin(2*pi*t)); % torque
tsim = 0:0.01:10; % simulation time

% ----- state space representation -----
% x1: theta
% x2: theta_dot
% x3: s
% x4: s_dot
f2 = @(t,X)((T(t) - c*X(2) + 2*k*l*X(3) + mp*g*l*sin(X(1)) - mp*l*X(3)*X(2)^2 - mp*g*(l*sin(X(1)) + X(3)*cos(X(1))) - mb*g*lc*sin(X(1)) - 2*mp*X(3)*X(4)*X(2))/...
    (mb*(rg^2+lc^2) + mp*(X(3)^2)) );
f4 = @(t,X)((-2*k*X(3) - mp*g*sin(X(1)) + mp*X(3)*X(2)^2 - mp*l*f2(t,X))/mp);

f = @(t,X)[X(2) ; f2(t,X) ; X(4) ; f4(t,X)];
% --------------------------------------

% ----- solve state equation -----
[ts,ys] = ode45(f,tsim,ic);
% --------------------------------

%% plot
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot(ts(:,1),ys(:,1));
title('theta-t plot');
xlabel('t(s)');
ylabel('theta(rad)');

subplot(2,2,2);
plot(ts(:,1),ys(:,2));
title('w-t plot');
xlabel('t(s)');
ylabel('w');

subplot(2,2,3);
plot(ts(:,1),ys(:,3));
title('s-t plot');
xlabel('t(s)');
ylabel('s');

subplot(2,2,4);
plot(ts(:,1),ys(:,4));
title('sdot-t plot');
xlabel('t(s)');
ylabel('sdot');

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
plot(ys(:,1),ys(:,3));
title({['phase portrait'];['theta - s']});
xlabel('s');
ylabel('theta');

subplot(2,1,2);
plot(ys(:,2),ys(:,4));
title({['phase portrait'];['w - sdot']});
xlabel('sdot');
ylabel('thetadot');