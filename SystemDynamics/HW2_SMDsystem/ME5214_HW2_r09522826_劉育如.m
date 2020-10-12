clear all; close all; clc;
% Homework2 matlab code
% main function of SMD system

%% Parameters
m1 = 1; c1 = 5; k1 = 6; % question1
m2 = 1; c2 = 2; k2 = 2; % question2
m3 = 1; c3 = 5; k3 = 6; % question3

% ----- sin function parameters -----
As = 10;
ws = 10;
% -----------------------------------

% ----- PWM function parameters -----
Ap = 10; 
alpha = 0.5;
wp = 10;
tp = 2*pi/wp;
wn = (k3/m3)^(1/2); % natural frequency
% -----------------------------------

t = 10; %sec

%% ODE45
% ----- state space representation -----
f1 = @(t,X)[X(2);1/m1*Force1(t,ws,As)-k1/m1*X(1)-c1/m1*X(2)]; % sin input m=1,c=5,k=6
f2 = @(t,X)[X(2);1/m2*Force1(t,ws,As)-k2/m2*X(1)-c2/m2*X(2)]; % sin input m=1,c=2,k=2
f3 = @(t,X)[X(2);1/m3*Force2(t,wp,alpha,Ap)-k3/m3*X(1)-c3/m3*X(2)]; % PWM input
% --------------------------------------

%% Question #1
[ts1,xs1] = ode45(f1,[0,t],[0;0]);

%% Question #2
[ts2,xs2] = ode45(f2,[0,t],[0;0]);

%% Question #3
[ts3,xs3] = ode45(f3,[0,t],[0;0]);

%% plot
figure('Name','Question #1','NumberTitle','off');
plot(ts1(:,1),xs1(:,1));
title('x-t plot');
xlabel('t');
ylabel('x');
txt = {'SMD Parameters:',['m=',num2str(m1),', c=',num2str(c1),', k=',num2str(k1)],'input: single sinusoidal wave',['->A=',num2str(As),', w=',num2str(ws)]};
text(max(ts1(:,1))-5,max(xs1(:,1))*0.8,txt)

figure('Name','Question #2','NumberTitle','off');
plot(ts2(:,1),xs2(:,1));
title('x-t plot');
xlabel('t');
ylabel('x');
txt = {'SMD Parameters:',['m=',num2str(m2),', c=',num2str(c2),', k=',num2str(k2)],'input: single sinusoidal wave',['->A=',num2str(As),', w=',num2str(ws)]};
text(max(ts1(:,1))-5,max(xs1(:,1))*0.8,txt)

figure('Name','Question #3','NumberTitle','off');
plot(ts3(:,1),xs3(:,1));
title('x-t plot');
xlabel('t');
ylabel('x');
txt = {'SMD Parameters:',['m=',num2str(m3),', c=',num2str(c3),', k=',num2str(k3)],'input: periodic squared wave',['->A=',num2str(Ap),' ,w=',num2str(wp),' ,alpha=',num2str(alpha)]};
text(max(ts1(:,1))-5,max(xs1(:,1))*0.8,txt)

%% Question #4
figure('Name','Question #4','NumberTitle','off');

% ----- Control Group -----
fc = @(t,X)[X(2);1/m3*Ap/2-k3/m3*X(1)-c3/m3*X(2)];
[ts_c,xs_c] = ode45(fc,[0,t],[0;0]);

subplot(2,1,1);
plot(ts_c(:,1),xs_c(:,1),'DisplayName','alpha = 1'); hold on;
title('x-t plot');
xlabel('t');
ylabel('x');
legend;

subplot(2,1,2);
plot(ts_c(:,1),xs_c(:,2),'DisplayName','alpha = 1'); hold on;
title('v-t plot')
xlabel('t');
ylabel('v');
legend;

% ----- Experimental Group ------
w = [wn wn*10 wn*100];
for i = 1:length(w)
    w3 = w(i);
    A = Ap;
    
    % ----- state space representation -----
    f3 = @(t,X)[X(2);1/m3*Force2(t,w3,alpha,A)-k3/m3*X(1)-c3/m3*X(2)]; % PWM input
    % --------------------------------------
    [ts,xs] = ode45(f3,[0,t],[0;0]);
    
    subplot(2,1,1);
    plot(ts(:,1),xs(:,1),'DisplayName',['w=',num2str(w3)]); hold on;
    title('x-t plot');
    xlabel('t');
    ylabel('x');
    legend;
    
    subplot(2,1,2);
    plot(ts(:,1),xs(:,2),'DisplayName',['w=',num2str(w3)]); hold on;
    title('v-t plot')
    xlabel('t');
    ylabel('v');
    legend;
end


%% one pulse sine function
function f = Force1(t,w,A)
    if t <= 2*pi/w;
        f = A*sin(w*t);
    else
        f = 0;
    end
end

%% PWM function
function f = Force2(t,w,alpha,A2)
    tf = mod(t,2*pi/w);
    if tf <= alpha*2*pi/w;
        f = A2;
    else
        f = 0;
    end
end