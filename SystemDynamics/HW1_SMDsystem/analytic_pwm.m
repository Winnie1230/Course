clear all; close all; clc;
% main function of SMD system
% compare different conditions(different w)
m = 1; c=5; k = 6;
t = 10; %sec

% ----- PWM function parameters -----
A2 = 10; 
alpha = 0.5;
w2 = 10;
tp = 2*pi/w2;
% -----------------------------------

% ----- state space representation -----
f2 = @(t,X)[X(2);1/m*Force2(t,w2,alpha,A2)-k/m*X(1)-c/m*X(2)]; % PWM input
% --------------------------------------
[ts,xs] = ode45(f2,[0,t],[0;0]);

ta = linspace(0,t,10000)';
%% analytic solution
y(:,1) = analytic(ta,tp,alpha,A2);

%% plot
plot(ts(:,1),xs(:,1),'DisplayName','ode45'); hold on;
plot(ta,y(:,1),'DisplayName','analytic solution'); hold on;
title('x-t plot');
% xlim([0 0.3]);
% ylim([-10 0]);
xlabel('t');
ylabel('x');
legend;

%% analytic function
function ant = analytic(t,tp,alpha,A)
    n=0;
    ant = zeros(length(t),1);
    for i = 1:3
        ant = ant + (1/6-1/2.*exp(-2.*(t-n.*tp))+1/3.*exp(-3.*(t-n.*tp))).*heaviside(t-n.*tp)-(1/6-1/2.*exp(-2.*(t-(alpha+n).*tp))+1/3.*exp(-3.*(t-(alpha+n).*tp))).*heaviside(t-(alpha+n).*tp);
        n = n + 1;
    end
    ant = A .* ant;
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