clear all; close all; clc;
% main function of SMD system
% compare different conditions(different w)
m = 1; c=5; k = 6;
t = 10; %sec
wn = (k/m)^(1/2); % natural frequency

% ----- PWM function parameters -----
A2 = 10; 
alpha2 = [0.2 0.5 0.7];
w2 = 10;
% -----------------------------------

% ----- Control Group -----
% fc = @(t,X)[X(2);1/m*A2/2-k/m*X(1)-c/m*X(2)];
% [ts_c,xs_c] = ode45(fc,[0,t],[0;0]);
% 
% subplot(2,1,1);
% plot(ts_c(:,1),xs_c(:,1),'DisplayName','alpha = 1'); hold on;
% title('x-t plot');
% xlabel('t');
% ylabel('x');
% legend;
% 
% subplot(2,1,2);
% plot(ts_c(:,1),xs_c(:,2),'DisplayName','alpha = 1'); hold on;
% title('v-t plot')
% xlabel('t');
% ylabel('v');
% legend;

% ----- Experimental Group ------
for i = 1:length(alpha2)
    alpha = alpha2(i);
    
    % ----- state space representation -----
    f2 = @(t,X)[X(2);1/m*Force2(t,w2,alpha,A2)-k/m*X(1)-c/m*X(2)]; % PWM input
    % --------------------------------------
    [ts,xs] = ode45(f2,[0,t],[0;0]);
    
    subplot(2,1,1);
    plot(ts(:,1),xs(:,1),'DisplayName',['alpha=',num2str(alpha)]); hold on;
    title('x-t plot');
    xlabel('t');
    ylabel('x');
    legend;
    
    subplot(2,1,2);
    plot(ts(:,1),xs(:,2),'DisplayName',['alpha=',num2str(alpha)]); hold on;
    title('v-t plot')
    xlabel('t');
    ylabel('v');
    legend;
end

%% PWM function
function f = Force2(t,w,alpha,A2)
    tf = mod(t,2*pi/w);
    if tf <= alpha*2*pi/w
        f = A2;
    else
        f = 0;
    end
end