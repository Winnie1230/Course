clear all; close all; clc;
% main function of SMD system animation
m = 1; c=5; k = 6;

% ----- sin function parameters -----
A1 = 20;
w1 = 5;
% -----------------------------------

equ_pt = [6,2.5];
t = 5; %sec

% ----- state space representation -----
f1 = @(t,X)[X(2);1/m*Force1(t,w1,A1)-k/m*X(1)-c/m*X(2)]; % sin input
% --------------------------------------

ic = [4;2];
[ts,xs] = ode45(f1,[0,t],ic);

%% animation
subplot(2,2,2);
h1 = animatedline('Color','b','LineWidth',1);
axis([min(ts(:,1)) max(ts(:,1)) -A1 A1]);
title('f-t');
xlabel('x');
ylabel('f');

subplot(2,2,3);
h2 = animatedline('Color','b','LineWidth',1);
axis([0 t min(xs(:,1)) max(xs(:,1))]);
title('x-t plot');
xlabel('t');
ylabel('x');

subplot(2,2,4);
h3 = animatedline('Color','b','LineWidth',1);
axis([0 t min(xs(:,2)) max(xs(:,2))]);
title('v-t plot');
xlabel('t');
ylabel('v');

%----- one pulse sine function input ------
for i=1:size(ts,1)
    subplot(2,2,1);
    PlotSMD([equ_pt(1)+xs(i),equ_pt(2)],equ_pt);

    addpoints(h1,ts(i,1),Force1(ts(i,1),w1,A1));
    addpoints(h2,ts(i,1),xs(i,1));
    addpoints(h3,ts(i,1),xs(i,2));
    drawnow;
end

%% one pulse sine function
function f = Force1(t,w,A)
    if t <= 2*pi/w;
        f = A*sin(w*t);
    else
        f = 0;
    end
end