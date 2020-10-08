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

%% analytic solution
c1 = -5*A1*w1/((w1^2+4)*(w1^2+9));
c2 = (-A1*w1^3+6*A1*w1)/((w1^2+4)*(w1^2+9));
c3 = A1*w1/(w1^2+4);
c4 = -A1*w1/(w1^2+9);

ta = linspace(0,t,1000)';
ant = c1.*cos(w1.*ta) + c2/w1.*sin(w1.*ta) + c3.*exp(-2.*ta) + c4.*exp(-3.*ta) - (c1.*cos(w1.*(ta-2.*pi/w1)) + c2/w1.*sin(w1.*(ta-2*pi/w1)) + c3.*exp(-2.*(ta-2*pi/w1)) + c4.*exp(-3.*(ta-2*pi/w1))).*heaviside(ta-2*pi/w1);

ant1 = @(t) c1*cos(w1*t) + c2/w1*sin(w1*t) + c3*exp(-2*t) + c4*exp(-3*t);
diff_ant1 = @(t) -w1*c1*sin(w1*t) + c2*cos(w1*t) -2*c3*exp(-2*t) -3*c4*exp(-3*t);

find = 0;
ic2 = [ant1(2*pi/w1),diff_ant1(2*pi/w1)];

c5 = -2*ic2(1)+(ic2(2)+5*ic2(1));
c6 = 3*ic2(1)-(ic2(2)+5*ic2(1));
ant2 = @(t) c5*exp(-2*t) + c6*exp(-3*t);

for i = 1:length(ta)
    if ta(i,1) <= 2*pi/w1
        y1(i,1) = ant1(ta(i,1));
    else
        y1(i,1) = ant2(ta(i,1)-2*pi/w1);
    end
end

%% plot
plot(ts(:,1),xs(:,1),'DisplayName','ode45'); hold on;
plot(ta(:,1),ant(:,1),'DisplayName','analytic solution - method 1'); hold on;
plot(ta(:,1),y1(:,1),'DisplayName','analytic solution - method 2'); hold on;
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