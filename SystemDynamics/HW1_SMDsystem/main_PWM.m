clear all; close all; clc;
% main function of SMD system
m = 1; c=5; k = 6;

% ----- PWM function parameters -----
A2 = 10; 
alpha = 0.6;
w2 = 5;
% -----------------------------------

equ_pt = [5,2.5];
t = 10; %sec

% ----- state space representation -----
f2 = @(t,X)[X(2);Force2(t,w2,alpha,A2)-6*X(1)-5*X(2)]; % PWM input
% --------------------------------------

[ts,xs] = ode45(f2,[0,t],[0;0]);

subplot(2,2,2);
h1 = animatedline('Color','b','LineWidth',1);
axis([0 t 0 A2]);
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

% ------ PWM input ------
for i=1:size(ts,1)
    subplot(2,2,1);
    PlotSMD([equ_pt(1)+xs(i),equ_pt(2)],equ_pt);

    addpoints(h1,ts(i,1),Force2(ts(i,1),w2,alpha,A2));
    addpoints(h2,ts(i,1),xs(i,1));
    addpoints(h3,ts(i,1),xs(i,2));
    drawnow;
end

% ------ PWM function ------
function f = Force2(t,w,alpha,A2)
    tf = mod(t,2*pi/w);
    if tf <= alpha*2*pi/w;
        f = A2;
    else
        f = 0;
    end
end