clear all; close all; clc;
% main function of SMD system
% analytic solution of m = 1, c = 5, k = 6;
% plot analytic solution and ode45 solution
m = 1; c = 5; k = 6;
t = 10; %sec

% ----- sin function parameters -----
A1 = 10;
w1 = 10;
% -----------------------------------

%% ODE45
% ----- state space representation -----
f1 = @(t,X)[X(2);1/m*Force1(t,w1,A1)-k/m*X(1)-c/m*X(2)]; % sin input
% --------------------------------------

[ts,xs] = ode45(f1,[0,t],[0;0]);

% ----- calculate force -----
for j = 1:size(ts,1)
    fs(j,1) = Force1(ts(j,1),w1,A1);
end

%% analytic solution
c1 = -5*A1*w1/((w1^2+4)*(w1^2+9));
c2 = (-A1*w1^3+6*A1*w1)/((w1^2+4)*(w1^2+9));
c3 = A1*w1/(w1^2+4);
c4 = -A1*w1/(w1^2+9);

ant1 = @(t) c1*cos(w1*t) + c2/w1*sin(w1*t) + c3*exp(-2*t) + c4*exp(-3*t);

find = 0;
ic2 = [0,0];
for i = 1:length(ts)
    if ts(i,1) > 2*pi/w1 & find == 0
        ic2 = [xs(i,1),xs(i,2)]; % initial condition of second part analytic solution
        break;
    end
end

c5 = -2*ic2(1)+(ic2(2)+5*ic2(1));
c6 = 3*ic2(1)-(ic2(2)+5*ic2(1));
ant2 = @(t) c5*exp(-2*t) + c6*exp(-3*t);

for i = 1:length(ts)
    if ts(i,1) <= 2*pi/w1
        y1(i,1) = ant1(ts(i,1));
    else
        y1(i,1) = ant2(ts(i,1)-2*pi/w1);
    end
end

%% plot
plot(ts(:,1),xs(:,1),'DisplayName','ode45'); hold on;
plot(ts(:,1),y1(:,1),'DisplayName','analytic solution'); hold on;
title('x-t plot');
xlabel('t');
ylabel('x');
legend;

%% one pulse sine function
function f = Force1(t,w,A)
    if t <= 2*pi/w;
        f = A*sin(w*t);
    else
        f = 0;
    end
end