clear all; close all; clc;
% main function of SMD system
% compare different conditions(different w and different A)
m = 1; c = 2; k = 2;
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

%% analytic solution
c1=-2*A1*w1/(w1^4+4);
c2 = -A1*w1*(w1^2-2)/(w1^4+4);
c3 = 2*A1*w1/(w1^4+4);
c4 = (A1*w1^4-4*A1*w1+4*A1)/(w1*(w1^4+4));

ant1 = @(t) c1*cos(w1*t) + c2/w1*sin(w1*t) + c3*exp(-t)*cos(t)+ c4*exp(-t)*sin(t);
% ant1 = @(t) c1*cos(w1*t) + c2/w1*sin(w1*t) + c3*exp(-t)*cos(w1*t) + c4*exp(-t)*sin(w1*t);

find = 0;
ic2 = [0,0];
for i = 1:length(ts)
    if ts(i,1) > 2*pi/w1 & find == 0
        ic2 = [xs(i,1),xs(i,2)]; % initial condition of second part analytic solution
        break;
    end
end

c5 = ic2(1);
c6 = ic2(1)+ic2(2);
ant2 = @(t) c5*exp(-t)*cos(t) + c6*exp(-t)*sin(t);

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